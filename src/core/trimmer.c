/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "amplicon.h"
#include "chimera.h"
#include "counter.h"
#include "fasta.h"
#include "fastq-reader.h"
#include "fastq-writer.h"
#include "lengths.h"
#include "logger.h"
#include "read-cover.h"
#include "reference.h"
#include "sam-tags.h"
#include "scheme.h"
#include "strbuf.h"
#include "test.h"
#include "trimmer.h"

struct trimmer_t {
  const trimmer_config_t *config;
  amplicon_analyzer_t *analyzer;
  chimera_finder_t *chimera;  /* NULL: the reads are not split */
  reference_set_t *refs;
  adapters_t *adapters;
  strbuf_t title;
  /* the cleavage site of the HA reference of the sample being trimmed
     (NULL: not read), set by trimmer_restrict_references() */
  const cleavage_reference_t *cleavage;
  /* the anchors of the markers of the sample being trimmed, one per
     marker (NULL where the sample has none), set by
     trimmer_restrict_references() */
  const marker_reference_t *markers[MARKER_N_MARKERS];
};

trimmer_t* trimmer_new(const trimmer_config_t *config) {
  assert(config != NULL);
  assert(config->scheme != NULL);
  trimmer_t *self;
  XMALLOC(self, sizeof(trimmer_t));
  self->config = config;
  self->analyzer = amplicon_analyzer_new(config->scheme, config->max_errors,
                                         config->policy, config->window);
  self->chimera = config->chimeras != READS_KEEP ?
      chimera_finder_new(config->scheme, config->max_errors,
                         config->window) : NULL;
  /* every reference allowed until a sample restricts them */
  self->cleavage = cleavage_anchors_chosen(config->cleavage, NULL);
  marker_anchors_chosen(config->markers, NULL, self->markers);
  self->refs = NULL;
  if (config->references != NULL) {
    assert(config->index != NULL);
    self->refs = reference_set_new(config->index, config->probe_len,
                                   config->probe_error_rate);
    reference_set_set_exhaustive(self->refs, config->exhaustive);
    reference_set_set_seeds_only(self->refs, config->seeds_only);
    reference_set_restrict_segment(self->refs, config->segment);
    amplicon_analyzer_use_references(self->analyzer, self->refs,
                                     config->terminal_window);
    if (self->chimera != NULL)
      chimera_finder_use_references(self->chimera, self->refs);
  }
  self->adapters = NULL;
  if (config->adapters != NULL) {
    self->adapters = adapters_new(config->adapters, config->window, NULL);
    amplicon_analyzer_use_adapters(self->analyzer, self->adapters);
  }
  strbuf_init(&self->title);
  return self;
}

void trimmer_restrict_virus(trimmer_t *self, virus_set_t types) {
  amplicon_analyzer_restrict_virus(self->analyzer, types);
  if (self->chimera != NULL)
    chimera_finder_restrict_virus(self->chimera, types);
}

void trimmer_restrict_references(trimmer_t *self, const bool *allowed) {
  amplicon_analyzer_restrict_references(self->analyzer, allowed);
  self->cleavage = cleavage_anchors_chosen(self->config->cleavage, allowed);
  marker_anchors_chosen(self->config->markers, allowed, self->markers);
}

void trimmer_delete(trimmer_t *self) {
  if (self == NULL)
    return;
  amplicon_analyzer_delete(self->analyzer);
  chimera_finder_delete(self->chimera);
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

/* append the title of a read to <t>, its name (the first word)
   followed by _<part> (<part> > 0: a part of a split read) and by
   _j<flank> (<flank> > 0: a flank of a read cut at its junction) */
static void append_title(strbuf_t *t, const char *title, size_t title_len,
    size_t part, size_t flank) {
  size_t name_len = 0;
  while (name_len < title_len && title[name_len] != ' ' &&
         title[name_len] != '\t')
    name_len++;
  strbuf_append(t, title, name_len);
  if (part > 0) {
    strbuf_append_char(t, '_');
    strbuf_append_long(t, (long)part);
  }
  if (flank > 0) {
    strbuf_append_str(t, "_j");
    strbuf_append_long(t, (long)flank);
  }
  strbuf_append(t, title + name_len, title_len - name_len);
}

/* the title of the record written for a read, or for the part <part>
   (from 1) of the <n_parts> parts of a split read, of length <len>: its
   title, then the tags; <category>: of the read (READ_CATEGORY_*);
   <flank>: the piece (from 1) of the <n_flanks> pieces of a read cut at
   its junction (and at the short deletions next to it), else 0 */
static void make_title(trimmer_t *self, const char *title, size_t title_len,
    const amplicon_t *a, size_t category, size_t len, size_t part,
    size_t n_parts, size_t flank, size_t n_flanks) {
  strbuf_t *t = &self->title;
  strbuf_clear(t);
  append_title(t, title, title_len, n_parts > 1 ? part : 0, flank);
  if (!self->config->tags)
    return;
  const char *strand = a->strand == STRAND_PLUS ? "+" :
                       a->strand == STRAND_MINUS ? "-" : ".";
  sam_tag_append_string(t, "fs", strand, 1);
  append_end_tags(self, &a->left, 'l', self->config->scheme->left);
  append_end_tags(self, &a->right, 'r', self->config->scheme->right);
  if (self->refs != NULL) {
    const char *name = read_category_name(category);
    sam_tag_append_string(t, "ct", name, strlen(name));
  }
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
  if (n_parts > 1) {
    sam_tag_start(t, "sp", 'Z');
    strbuf_append_long(t, (long)part);
    strbuf_append_char(t, '/');
    strbuf_append_long(t, (long)n_parts);
  }
  if (flank > 0) {
    sam_tag_start(t, "jf", 'Z');
    strbuf_append_long(t, (long)flank);
    strbuf_append_char(t, '/');
    strbuf_append_long(t, (long)n_flanks);
  }
}

/* the bounds of the pieces of the kept part of the read <a> cut at its
   junction <cut> and at the short deletions next to it, in the order of
   the read: the first base of each piece, then the base after the last
   one, into <bounds>; their number */
static size_t piece_bounds(const amplicon_t *a, long cut, long *bounds) {
  /* the cut points in increasing order (few: an insertion sort) */
  long points[FLANK_DELETIONS_MAX_CUTS + 1];
  size_t n_points = 0;
  points[n_points++] = cut;
  for (size_t i = 0; i < a->n_flank_cuts; i++) {
    size_t j = n_points++;
    while (j > 0 && points[j - 1] > a->flank_cuts[i]) {
      points[j] = points[j - 1];
      j--;
    }
    points[j] = a->flank_cuts[i];
  }
  /* those inside the kept part, each once */
  size_t n = 0;
  bounds[n++] = a->keep_from;
  for (size_t i = 0; i < n_points; i++)
    if (points[i] > bounds[n - 1] && points[i] <= a->keep_to)
      bounds[n++] = points[i];
  bounds[n++] = a->keep_to + 1;
  return n;
}

/* analyze, count and write the bases <from> to <to> of the read <rec>
   of <batch>, the part <part> (from 1) of <n_parts>; the record is
   appended to the buffer of <out> of its category, or to the only one
   (see trimmer_n_outputs()) */
/*
   Read the cleavage site of the haemagglutinin out of a read
   (cleavage-site.h). Which reads are read differs between the
   platforms, because what reaches the site differs:

   - reads of whole amplicons (Nanopore) reach the site only when they
     are full templates, as the internally primed ones begin just after
     it and the genomes with an internal deletion have removed it; they
     are located on a reference, so only those on the HA reference of
     the sample are searched;
   - a read of a fragment (Illumina) which covers the site carries no
     primer, so it has no category and is not located at all: every read
     is searched, and the anchor itself says which ones cover the site.
*/
static void read_cleavage_site(trimmer_t *self, const char *seq, size_t len,
    const read_cover_t *cover, size_t category, counter_t *counts) {
  if (self->cleavage == NULL)
    return;
  if (!self->config->illumina &&
      !(cover->located && cover->ref == self->cleavage->anchor.ref &&
        category == READ_CATEGORY_FULL))
    return;
  counter_set_cleavage_reference(counts, self->cleavage->aa,
                                 self->cleavage->fusion);
  counter_add_cleavage_searched(counts);
  /* one search of the anchor per read; the site and the control
     windows are then read from either side of the match */
  anchor_hit_t hit;
  if (!anchor_find(&self->cleavage->anchor, seq, len,
                   self->config->probe_error_rate, &hit))
    return;
  char window[CLEAVAGE_SITE_NT + 1];
  if (anchor_bases_before(&hit, seq, len, CLEAVAGE_SITE_NT, window))
    counter_add_cleavage(counts, window);
  if (anchor_bases_at(&hit, seq, len, CLEAVAGE_ANCHOR_NT, window))
    counter_add_cleavage_fusion(counts, window);
  for (size_t i = 0; i < CLEAVAGE_N_CONTROLS; i++)
    if (anchor_bases_after(&hit, seq, len, i * CLEAVAGE_SITE_NT,
                           CLEAVAGE_SITE_NT, window))
      counter_add_cleavage_control(counts, i, window);
}

/* analyze, count and write a part of a read (the whole read, or one of
   its amplicons); <write>: false for a part of a read with several
   amplicons filtered out, which is counted as such */
static void trim_part(trimmer_t *self, const fastq_batch_t *batch,
    const fastq_record_t *rec, long from, long to, size_t part,
    size_t n_parts, bool write, strbuf_t *out, counter_t *counts) {
  const char *data = batch->data.s;
  const char *seq = data + rec->seq_from + from;
  size_t len = (size_t)(to - from + 1);
  const amplicon_t *a = amplicon_analyzer_run(self->analyzer, seq, len);
  read_cover_t cover;
  read_cover_locate(self->config->references, self->refs, a, &cover);
  size_t category = read_cover_category(&cover);
  counter_add_located_read(counts, a, len, &cover);
  read_cleavage_site(self, seq, len, &cover, category, counts);
  /* the markers of the segment of the read; a read which is not
     located (the reads of a fragment, see marker-read.h) is searched
     for every marker. A read gives only the markers which lie whole on
     a stretch where it follows its reference: between the cuts of an
     internal deletion, and clear of a primer which matched inside it,
     whose bases are the primer and not the template. A read with an
     internal deletion whose junction was not located gives none, as
     which of its bases belong together is then not known */
  long cuts[1 + FLANK_DELETIONS_MAX_CUTS];
  marker_blocked_t blocked[AMPLICON_MAX_INTERNAL];
  marker_breaks_t breaks = {cuts, 0, blocked, 0};
  bool readable = category != READ_CATEGORY_DELETION || a->junction.found;
  if (readable) {
    if (category == READ_CATEGORY_DELETION) {
      cuts[breaks.n_cuts++] = a->junction_at;
      for (size_t i = 0; i < a->n_flank_cuts; i++)
        cuts[breaks.n_cuts++] = a->flank_cuts[i];
    }
    for (size_t i = 0; i < a->n_internal; i++) {
      blocked[breaks.n_blocked].from = a->internal[i].start;
      blocked[breaks.n_blocked].to = a->internal[i].end;
      breaks.n_blocked++;
    }
    counter_add_markers(counts, self->markers,
                        cover.located ? (long)cover.ref : -1, &breaks,
                        seq, len, self->config->probe_error_rate);
  }
  if (!write || (self->config->filtered != NULL &&
                 self->config->filtered[category])) {
    counter_add_filtered(counts, category);
    return;
  }
  size_t kept = amplicon_kept_length(a);
  if (kept == 0 || kept < self->config->min_length) {
    counter_add_too_short(counts);
    return;
  }
  if (self->config->by_category)
    out += category;
  const char *qual = data + rec->qual_from + from;
  long cut = self->config->deletions == READS_CUT &&
      category == READ_CATEGORY_DELETION ? amplicon_junction_cut(a) : -1;
  if (cut < 0) {
    make_title(self, data + rec->title_from, rec->title_len, a, category,
               len, part, n_parts, 0, 0);
    fastq_writer_record(out, self->title.s, self->title.len,
                        seq + a->keep_from, qual + a->keep_from, kept);
    counter_add_written(counts, kept, category, rec->qual_missing ? NULL :
                        qual + a->keep_from);
    return;
  }
  /* the pieces of the kept part, in the order of the read: its two
     flanks, cut at the junction, and cut again at the short deletions
     next to it */
  counter_add_cut(counts);
  long bounds[FLANK_DELETIONS_MAX_CUTS + 3];
  size_t n_flanks = piece_bounds(a, cut, bounds) - 1;
  if (n_flanks > 2)
    counter_add_cut_again(counts);
  for (size_t flank = 1; flank <= n_flanks; flank++) {
    long flank_from = bounds[flank - 1];
    size_t flank_len = (size_t)(bounds[flank] - flank_from);
    if (flank_len < self->config->min_length) {
      counter_add_too_short(counts);
      continue;
    }
    make_title(self, data + rec->title_from, rec->title_len, a, category,
               len, part, n_parts, flank, n_flanks);
    fastq_writer_record(out, self->title.s, self->title.len,
                        seq + flank_from, qual + flank_from, flank_len);
    counter_add_written(counts, flank_len, category,
                        rec->qual_missing ? NULL : qual + flank_from);
  }
}

void trimmer_run(trimmer_t *self, const fastq_batch_t *batch, size_t from,
    size_t to, strbuf_t *out, counter_t *counts) {
  assert(to <= batch->n_records);
  const char *data = batch->data.s;
  for (size_t i = from; i < to; i++) {
    const fastq_record_t *rec = batch->records + i;
    /* (the title of the first read of the range is read for the metadata
       of the run) */
    counter_add_input(counts, rec->seq_len, rec->qual_missing ? NULL :
                      data + rec->qual_from,
                      i == from ? data + rec->title_from : NULL,
                      rec->title_len);
    chimera_part_t parts[CHIMERA_MAX_PARTS];
    size_t n_parts = self->chimera == NULL ? 0 :
        chimera_finder_run(self->chimera, data + rec->seq_from,
                           rec->seq_len, parts);
    if (n_parts < 2) {
      /* (an empty read is analyzed as well, and counted) */
      trim_part(self, batch, rec, 0, (long)rec->seq_len - 1, 1, 1, true,
                out, counts);
      continue;
    }
    size_t parts_len = 0;
    bool write = self->config->chimeras != READS_FILTER;
    for (size_t p = 0; p < n_parts; p++) {
      trim_part(self, batch, rec, parts[p].from, parts[p].to, p + 1,
                n_parts, write, out, counts);
      parts_len += (size_t)(parts[p].to - parts[p].from + 1);
    }
    counter_add_split(counts, n_parts, rec->seq_len, parts_len);
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

/* the category of the reads: the tag, and the records written apart by
   category */
static void test_categories(const scheme_t *sc) {
  /* the interior of a segment (not similar to the primers) */
  const char *interior =
    "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG"
    "CTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGAC"
    "TGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGACAGGTCACGCAGAGGC"
    "GCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAA"
    "CTCCAGCGCGGTCAGTTCCATCACCCTAAGTAACCGAATAATGCGTTCGCTCTATTGACT";
  const char *uni12 = "AGCAAAAGCAGG";
  const char *rc13 = "CCTTGTTTCTACTGATCACGCGT";
  char segment[512], read[512];
  snprintf(segment, sizeof(segment), "%s%sCCTTGTTTCTACT", uni12, interior);
  fasta_record_t record = {.name = "seg", .description = "",
                           .seq = segment, .seq_len = strlen(segment)};
  fasta_t fa = {.records = &record, .n_records = 1};
  reference_index_t *index = reference_index_new(&fa);
  trimmer_config_t config = {.scheme = sc, .max_errors = 0.1,
      .policy = TRIM_TERMINUS, .window = AMPLICON_DEFAULT_WINDOW,
      .references = &fa, .index = index,
      .probe_len = REFERENCE_DEFAULT_PROBE_LEN,
      .probe_error_rate = REFERENCE_DEFAULT_MAX_ERROR_RATE,
      .terminal_window = REFERENCE_DEFAULT_TERMINAL_WINDOW,
      .min_length = 1, .tags = true};
  test_eq_size(trimmer_n_outputs(&config), 1);
  trimmer_t *tr = trimmer_new(&config);
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  strbuf_t out[TRIMMER_MAX_OUTPUTS];
  for (size_t k = 0; k < TRIMMER_MAX_OUTPUTS; k++)
    strbuf_init(out + k);
  counter_t *counts = counter_new(sc, &fa);
  /* a full template, then an internally primed one: both in the only
     output, each with the tag of its category */
  snprintf(read, sizeof(read), "%s%s%s", uni12, interior, rc13);
  one_record(&batch, "whole", read);
  trimmer_run(tr, &batch, 0, 1, out, counts);
  snprintf(read, sizeof(read), "%s%s%s", uni12, interior + 150, rc13);
  one_record(&batch, "inside", read);
  trimmer_run(tr, &batch, 0, 1, out, counts);
  strbuf_append_char(out, '\0');
  test_assert(strstr(out[0].s, "\tct:Z:full\ttl:i:") != NULL);
  test_assert(strstr(out[0].s, "\tct:Z:internal\ttl:i:") != NULL);
  test_eq_size(out[1].len, 0);
  trimmer_delete(tr);
  /* written apart: each read in the output of its category */
  config.by_category = true;
  test_eq_size(trimmer_n_outputs(&config), READ_N_CATEGORIES);
  tr = trimmer_new(&config);
  strbuf_clear(out);
  trimmer_run(tr, &batch, 0, 1, out, counts);
  snprintf(read, sizeof(read), "%s%s%s", uni12, interior, rc13);
  one_record(&batch, "whole", read);
  trimmer_run(tr, &batch, 0, 1, out, counts);
  test_eq_mem(out[READ_CATEGORY_FULL].s, 7, "@whole\t");
  test_eq_mem(out[READ_CATEGORY_INTERNAL].s, 8, "@inside\t");
  test_eq_size(out[READ_CATEGORY_DELETION].len, 0);
  test_eq_size(out[READ_CATEGORY_FRAGMENT].len, 0);
  test_eq_size(out[READ_CATEGORY_UNKNOWN].len, 0);
  /* a read with an internal deletion (the bases 90 to 209 of the
     interior are missing): written whole */
  char flanks[512];
  snprintf(flanks, sizeof(flanks), "%.90s%s", interior, interior + 210);
  snprintf(read, sizeof(read), "%s%s%s", uni12, flanks, rc13);
  one_record(&batch, "deleted x=1", read);
  trimmer_run(tr, &batch, 0, 1, out, counts);
  strbuf_t *del = out + READ_CATEGORY_DELETION;
  test_eq_mem(del->s, 13, "@deleted x=1\t");
  strbuf_append_char(del, '\0');
  test_assert(strstr(del->s, "\tjf:Z:") == NULL);
  test_assert(strstr(del->s, flanks) != NULL);
  test_eq_size(counter_reads(counts)->cut, 0);
  trimmer_delete(tr);
  /* cut at the junction: two records, the two flanks, which together
     are the kept part of the read, the template with the conserved
     termini (the junction can be a few bases from
     where it was made, if the bases beside it are the same on its two
     sides) */
  char kept[512];
  snprintf(kept, sizeof(kept), "%s%sCCTTGTTTCTACT", uni12, flanks);
  config.deletions = READS_CUT;
  tr = trimmer_new(&config);
  strbuf_clear(del);
  size_t written = counter_reads(counts)->written[READ_CATEGORY_DELETION];
  trimmer_run(tr, &batch, 0, 1, out, counts);
  strbuf_append_char(del, '\0');
  test_eq_mem(del->s, 16, "@deleted_j1 x=1\t");
  const char *second = strstr(del->s, "\n@deleted_j2 x=1\t");
  test_assert(second != NULL);
  test_assert(strstr(del->s, "\tct:Z:deletion\t") != NULL);
  test_assert(strstr(del->s, "\tjf:Z:1/2\n") != NULL);
  test_assert(second != NULL && strstr(second, "\tjf:Z:2/2\n") != NULL);
  if (second != NULL) {
    const char *seq1 = strchr(del->s, '\n') + 1,
               *seq2 = strchr(second + 1, '\n') + 1;
    size_t len1 = strcspn(seq1, "\n"), len2 = strcspn(seq2, "\n");
    char joined[512];
    snprintf(joined, sizeof(joined), "%.*s%.*s", (int)len1, seq1,
             (int)len2, seq2);
    test_eq_size(len1 + len2, strlen(kept));
    test_eq_mem(joined, strlen(kept), kept);
    test_assert(len1 >= 97 && len1 <= 107);
    /* the qualities are cut with the bases */
    test_eq_mem(seq1 + len1 + 1, 2, "+\n");
    test_eq_size(strcspn(seq1 + len1 + 3, "\n"), len1);
  }
  test_eq_size(counter_reads(counts)->cut, 1);
  test_eq_size(counter_reads(counts)->written[READ_CATEGORY_DELETION],
               written + 2);
  /* the reads of the other categories are written as without the
     option */
  snprintf(read, sizeof(read), "%s%s%s", uni12, interior, rc13);
  one_record(&batch, "whole", read);
  strbuf_clear(out + READ_CATEGORY_FULL);
  trimmer_run(tr, &batch, 0, 1, out, counts);
  test_eq_mem(out[READ_CATEGORY_FULL].s, 7, "@whole\t");
  test_eq_size(counter_reads(counts)->cut, 1);
  /* the read of the minus strand: the flank 1 is the one before the
     junction in the read, the reverse complement of the right flank */
  char minus[512];
  snprintf(read, sizeof(read), "%s%s%s", uni12, flanks, rc13);
  size_t read_len = strlen(read);
  for (size_t i = 0; i < read_len; i++) {
    char b = read[read_len - 1 - i];
    minus[i] = b == 'A' ? 'T' : b == 'C' ? 'G' : b == 'G' ? 'C' : 'A';
  }
  minus[read_len] = '\0';
  one_record(&batch, "deleted", minus);
  strbuf_clear(del);
  trimmer_run(tr, &batch, 0, 1, out, counts);
  strbuf_append_char(del, '\0');
  test_eq_mem(del->s, 12, "@deleted_j1\t");
  second = strstr(del->s, "\n@deleted_j2\t");
  test_assert(second != NULL);
  if (second != NULL) {
    const char *seq1 = strchr(del->s, '\n') + 1,
               *seq2 = strchr(second + 1, '\n') + 1;
    size_t len1 = strcspn(seq1, "\n"), len2 = strcspn(seq2, "\n");
    test_eq_size(len1 + len2, strlen(kept));
    test_assert(len1 >= 98 && len1 <= 108);
    /* (the flank 2 ends with the reverse complement of the first
       bases of the segment) */
    test_eq_mem(seq2 + len2 - 12, 12, "CCTGCTTTTGCT");
  }
  test_eq_size(counter_reads(counts)->cut, 2);
  test_eq_size(counter_reads(counts)->cut_again, 0);
  /* a short deletion next to the junction too (the bases 50 to 59 of the
     interior, 30 bases before it): three pieces, cut at both */
  char flanks2[512];
  snprintf(flanks2, sizeof(flanks2), "%.50s%.30s%s", interior,
           interior + 60, interior + 210);
  snprintf(read, sizeof(read), "%s%s%s", uni12, flanks2, rc13);
  one_record(&batch, "twice", read);
  strbuf_clear(del);
  trimmer_run(tr, &batch, 0, 1, out, counts);
  strbuf_append_char(del, '\0');
  test_eq_mem(del->s, 10, "@twice_j1\t");
  test_assert(strstr(del->s, "\tjf:Z:1/3\n") != NULL);
  test_assert(strstr(del->s, "\n@twice_j2\t") != NULL);
  test_assert(strstr(del->s, "\n@twice_j3\t") != NULL);
  test_assert(strstr(del->s, "\tjf:Z:3/3\n") != NULL);
  test_eq_size(counter_reads(counts)->cut, 3);
  test_eq_size(counter_reads(counts)->cut_again, 1);
  trimmer_delete(tr);
  /* filtered out: counted, not written, the other categories written */
  bool filtered[READ_N_CATEGORIES] = {false};
  filtered[READ_CATEGORY_DELETION] = true;
  config.deletions = READS_FILTER;
  config.filtered = filtered;
  tr = trimmer_new(&config);
  strbuf_clear(del);
  trimmer_run(tr, &batch, 0, 1, out, counts);
  test_eq_size(del->len, 0);
  test_eq_size(counter_reads(counts)->filtered[READ_CATEGORY_DELETION], 1);
  snprintf(read, sizeof(read), "%s%s%s", uni12, interior, rc13);
  one_record(&batch, "whole", read);
  strbuf_clear(out + READ_CATEGORY_FULL);
  trimmer_run(tr, &batch, 0, 1, out, counts);
  test_eq_mem(out[READ_CATEGORY_FULL].s, 7, "@whole\t");
  trimmer_delete(tr);
  counter_delete(counts);
  for (size_t k = 0; k < TRIMMER_MAX_OUTPUTS; k++)
    strbuf_free(out + k);
  fastq_batch_free(&batch);
  reference_index_delete(index);
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

  /* two amplicons in a read: one record for each, the name followed
     by the number of the part */
  config.chimeras = READS_CUT;
  trimmer_t *split = trimmer_new(&config);
  snprintf(read, sizeof(read), "%s%s%s%s%s%s", uni12, body, rc13, "ACGT",
           uni12, body);
  one_record(&batch, "r1b x=1", read);
  strbuf_clear(&out);
  trimmer_run(split, &batch, 0, 1, &out, counts);
  snprintf(expected, sizeof(expected),
           "@r1b_1 x=1\tfs:Z:+\tfl:Z:uni12r\tfr:Z:mbtuni-13\ttl:i:12"
           "\ttr:i:13\tsp:Z:1/2\n%s\n+\n", body);
  test_eq_mem(out.s, strlen(expected), expected);
  test_assert(strstr(out.s, "\n@r1b_2 x=1\tfs:Z:+\tfl:Z:uni12r\ttl:i:12"
                            "\ttr:i:0\tsp:Z:2/2\n") != NULL);
  test_eq_size(stats->reads_in, 2);
  test_eq_size(stats->reads_out, 3);
  test_eq_size(stats->split, 1);
  /* (the bases of the input: those between the parts included) */
  test_eq_size(stats->bases_in, strlen(uni12) + strlen(body) +
               strlen(rc13) + strlen(read));
  trimmer_delete(split);
  /* filtered out: the parts analyzed and counted, none written */
  config.chimeras = READS_FILTER;
  split = trimmer_new(&config);
  counter_t *filter_counts = counter_new(sc, NULL);
  strbuf_clear(&out);
  trimmer_run(split, &batch, 0, 1, &out, filter_counts);
  test_eq_size(out.len, 0);
  test_eq_size(counter_reads(filter_counts)->split, 1);
  test_eq_size(counter_reads(filter_counts)->filtered[READ_CATEGORY_UNKNOWN],
               2);
  test_eq_size(counter_reads(filter_counts)->reads_out, 0);
  counter_delete(filter_counts);
  trimmer_delete(split);
  config.chimeras = READS_KEEP;

  /* a primer dimer is not written */
  snprintf(read, sizeof(read), "%s%s", uni12, rc13);
  one_record(&batch, "r2", read);
  strbuf_clear(&out);
  trimmer_run(tr, &batch, 0, 1, &out, counts);
  test_eq_size(out.len, 0);
  test_eq_size(stats->dimers, 1);
  test_eq_size(stats->too_short, 1);
  test_eq_size(stats->reads_in, 3);
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
  test_categories(sc);
  scheme_delete(sc);
  delete_logger(logger);
}
