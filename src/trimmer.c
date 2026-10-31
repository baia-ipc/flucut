/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "amplicon.h"
#include "counter.h"
#include "fasta.h"
#include "fastq-reader.h"
#include "fastq-writer.h"
#include "lengths.h"
#include "logger.h"
#include "reference.h"
#include "sam-tags.h"
#include "scheme.h"
#include "strbuf.h"
#include "test.h"
#include "trimmer.h"

struct trimmer_t {
  const trimmer_config_t *config;
  amplicon_analyzer_t *analyzer;
  reference_set_t *refs;
  adapters_t *adapters;
  strbuf_t title;
};

trimmer_t* trimmer_new(const trimmer_config_t *config) {
  assert(config != NULL);
  assert(config->scheme != NULL);
  trimmer_t *self;
  XMALLOC(self, sizeof(trimmer_t));
  self->config = config;
  self->analyzer = amplicon_analyzer_new(config->scheme, config->max_errors,
                                         config->policy, config->window);
  self->refs = NULL;
  if (config->references != NULL) {
    assert(config->index != NULL);
    self->refs = reference_set_new(config->index, config->probe_len,
                                   config->probe_error_rate);
    reference_set_set_exhaustive(self->refs, config->exhaustive);
    reference_set_set_seeds_only(self->refs, config->seeds_only);
    amplicon_analyzer_use_references(self->analyzer, self->refs,
                                     config->terminal_window);
  }
  self->adapters = NULL;
  if (config->adapters != NULL) {
    self->adapters = adapters_new(config->adapters, config->window, NULL);
    amplicon_analyzer_use_adapters(self->analyzer, self->adapters);
  }
  strbuf_init(&self->title);
  return self;
}

void trimmer_restrict_subtype(trimmer_t *self, int subtype) {
  amplicon_analyzer_restrict_subtype(self->analyzer, subtype);
}

void trimmer_delete(trimmer_t *self) {
  if (self == NULL)
    return;
  amplicon_analyzer_delete(self->analyzer);
  reference_set_delete(self->refs);
  adapters_delete(self->adapters);
  strbuf_free(&self->title);
  free(self);
}

static const char *origin_names[] = {"unknown", "terminal", "internal"};

/* append the tags of an end (<l>: 'l' or 'r'; <primers>: the primers of
   that end of the scheme) */
static void append_end_tags(trimmer_t *self, const amplicon_end_t *end,
    char l, const primer_t *primers) {
  if (!end->found)
    return;
  strbuf_t *t = &self->title;
  char tag[3] = {'f', l, '\0'};
  const char *name = primers[end->variant].name;
  sam_tag_append_string(t, tag, name, strlen(name));
  if (self->refs == NULL)
    return;
  tag[0] = 'o';
  const char *origin = origin_names[end->origin];
  sam_tag_append_string(t, tag, origin, strlen(origin));
  if (!end->ref.mapped)
    return;
  tag[0] = 'm';
  sam_tag_start(t, tag, 'Z');
  strbuf_append_str(t,
      self->config->references->records[end->ref.ref_index].name);
  strbuf_append_char(t, ':');
  strbuf_append_long(t, end->ref.ref_start + 1);
  strbuf_append_char(t, ':');
  strbuf_append_char(t, end->ref.reverse ? '-' : '+');
}

/* the title of the record written for a read: its title, then the tags */
static void make_title(trimmer_t *self, const char *title, size_t title_len,
    const amplicon_t *a, size_t len) {
  strbuf_t *t = &self->title;
  strbuf_clear(t);
  strbuf_append(t, title, title_len);
  if (!self->config->tags)
    return;
  const char *strand = a->strand == STRAND_PLUS ? "+" :
                       a->strand == STRAND_MINUS ? "-" : ".";
  sam_tag_append_string(t, "fs", strand, 1);
  append_end_tags(self, &a->left, 'l', self->config->scheme->left);
  append_end_tags(self, &a->right, 'r', self->config->scheme->right);
  sam_tag_append_int(t, "tl", a->keep_from);
  sam_tag_append_int(t, "tr", (long)len - 1 - a->keep_to);
  if (a->n_internal > 0)
    sam_tag_append_int(t, "ni", (long)a->n_internal);
  if (a->adapter_start.found) {
    const char *kit = adapters_kit_name(a->adapter_start.kit);
    sam_tag_append_string(t, "a5", kit, strlen(kit));
  }
  if (a->adapter_end.found) {
    const char *kit = adapters_kit_name(a->adapter_end.kit);
    sam_tag_append_string(t, "a3", kit, strlen(kit));
  }
}

void trimmer_run(trimmer_t *self, const fastq_batch_t *batch, size_t from,
    size_t to, strbuf_t *out, counter_t *counts) {
  assert(to <= batch->n_records);
  const char *data = batch->data.s;
  for (size_t i = from; i < to; i++) {
    const fastq_record_t *rec = batch->records + i;
    const char *seq = data + rec->seq_from;
    const amplicon_t *a = amplicon_analyzer_run(self->analyzer, seq,
                                                rec->seq_len);
    counter_add_read(counts, a, rec->seq_len);
    size_t kept = amplicon_kept_length(a);
    if (kept == 0 || kept < self->config->min_length) {
      counter_add_too_short(counts);
      continue;
    }
    make_title(self, data + rec->title_from, rec->title_len, a,
               rec->seq_len);
    fastq_writer_record(out, self->title.s, self->title.len,
                        seq + a->keep_from,
                        data + rec->qual_from + a->keep_from, kept);
    counter_add_written(counts, kept);
  }
}

/* a batch of one record */
static void one_record(fastq_batch_t *batch, const char *title,
    const char *seq) {
  strbuf_clear(&batch->data);
  size_t len = strlen(seq);
  if (batch->records_cap == 0) {
    XMALLOC(batch->records, sizeof(fastq_record_t));
    batch->records_cap = 1;
  }
  fastq_record_t *r = batch->records;
  r->title_from = 0;
  r->title_len = strlen(title);
  strbuf_append_str(&batch->data, title);
  r->seq_from = batch->data.len;
  r->seq_len = len;
  strbuf_append_str(&batch->data, seq);
  r->qual_from = batch->data.len;
  for (size_t i = 0; i < len; i++)
    strbuf_append_char(&batch->data, (char)('!' + (char)(i % 40)));
  r->qual_missing = false;
  batch->n_records = 1;
}

void trimmer_test(void) {
  Logger *logger = new_logger(stderr, LOG_ERROR);
  scheme_t *sc = scheme_new_builtin(SCHEME_DEFAULT, logger);
  const char *uni12 = "AGCAAAAGCAGG";
  const char *rc13 = "CCTTGTTTCTACTGATCACGCGT";
  const char *body = "GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAA";
  trimmer_config_t config = {.scheme = sc, .max_errors = 0.1,
      .policy = TRIM_PRIMER, .window = AMPLICON_DEFAULT_WINDOW,
      .references = NULL, .min_length = 1, .tags = true};
  trimmer_t *tr = trimmer_new(&config);
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  strbuf_t out;
  strbuf_init(&out);
  counter_t *counts = counter_new(sc, NULL);
  const counter_reads_t *stats = counter_reads(counts);
  char read[256], expected[512];

  /* a read with both primers: the body is written, the qualities are
     cut with the sequence, the tags tell what was removed */
  snprintf(read, sizeof(read), "%s%s%s", uni12, body, rc13);
  one_record(&batch, "r1 x=1", read);
  trimmer_run(tr, &batch, 0, 1, &out, counts);
  snprintf(expected, sizeof(expected),
           "@r1 x=1\tfs:Z:+\tfl:Z:uni12r\tfr:Z:mbtuni-13\ttl:i:12"
           "\ttr:i:23\n%s\n+\n", body);
  test_eq_mem(out.s, strlen(expected), expected);
  test_eq_mem(out.s + strlen(expected), 3, "-./");
  test_eq_size(stats->reads_out, 1);
  test_eq_size(stats->bases_out, strlen(body));
  test_eq_size(stats->n_ends[2], 1);
  test_eq_size(lengths_count(counter_lengths(counts), 0,
                             strlen(body) / LENGTHS_BIN), 1);

  /* a primer dimer is not written */
  snprintf(read, sizeof(read), "%s%s", uni12, rc13);
  one_record(&batch, "r2", read);
  strbuf_clear(&out);
  trimmer_run(tr, &batch, 0, 1, &out, counts);
  test_eq_size(out.len, 0);
  test_eq_size(stats->dimers, 1);
  test_eq_size(stats->too_short, 1);
  test_eq_size(stats->reads_in, 2);
  trimmer_delete(tr);

  /* without tags, and with a minimum length */
  config.tags = false;
  config.min_length = strlen(body) + 1;
  tr = trimmer_new(&config);
  snprintf(read, sizeof(read), "%s%s%s", uni12, body, rc13);
  one_record(&batch, "r3", read);
  strbuf_clear(&out);
  trimmer_run(tr, &batch, 0, 1, &out, counts);
  test_eq_size(out.len, 0);
  config.min_length = strlen(body);
  trimmer_delete(tr);
  tr = trimmer_new(&config);
  trimmer_run(tr, &batch, 0, 1, &out, counts);
  snprintf(expected, sizeof(expected), "@r3\n%s\n+\n", body);
  test_eq_mem(out.s, strlen(expected), expected);
  trimmer_delete(tr);

  counter_delete(counts);
  strbuf_free(&out);
  fastq_batch_free(&batch);
  scheme_delete(sc);
  delete_logger(logger);
}
