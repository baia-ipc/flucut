/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdlib.h>
#include "alloc.h"
#include "amplicon.h"
#include "auto-settings.h"
#include "fastq-reader.h"
#include "numutils.h"
#include "reads-reader.h"
#include "reference.h"
#include "segments.h"
#include "strutils.h"
#include "test.h"

void auto_settings_init(auto_settings_t *self, const options_t *opts) {
  assert(self != NULL);
  assert(opts != NULL);
  bool automatic = opts->max_errors == OPTIONS_MAX_ERRORS_AUTO;
  self->max_errors = (report_setting_t){.automatic = automatic,
      .pending = automatic,
      .value = automatic ? AUTO_MAX_ERRORS_LOW : opts->max_errors,
      .measured = 0, .n_reads = 0, .file = NULL};
  self->file = NULL;
}

bool auto_settings_pending(const auto_settings_t *self) {
  assert(self != NULL);
  return self->max_errors.pending;
}

double auto_settings_max_errors(double error_rate) {
  return error_rate >= AUTO_MAX_ERRORS_HIGH_RATE ? AUTO_MAX_ERRORS_HIGH :
         AUTO_MAX_ERRORS_LOW;
}

/* add the errors and the aligned read bases of an end found */
static void add_end(const amplicon_end_t *end, size_t *errors,
    size_t *bases) {
  if (!end->found)
    return;
  *errors += (size_t)end->hit.n_errors;
  *bases += (size_t)(end->hit.end - end->hit.start + 1);
}

double auto_settings_error_rate(const trimmer_config_t *config,
    const char *file, size_t max_reads, size_t *n_reads, Logger *logger) {
  assert(config != NULL);
  assert(file != NULL);
  amplicon_analyzer_t *a = amplicon_analyzer_new(config->scheme,
      AUTO_MAX_ERRORS_LOW, config->policy, config->window);
  reads_reader_t *reader = reads_reader_new(file, logger);
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  size_t n = 0, errors = 0, bases = 0, got;
  while (n < max_reads &&
         (got = reads_reader_read_batch(reader, &batch, max_reads - n)) > 0) {
    for (size_t i = 0; i < got; i++) {
      const fastq_record_t *rec = batch.records + i;
      const amplicon_t *r = amplicon_analyzer_run(a,
          batch.data.s + rec->seq_from, rec->seq_len);
      add_end(&r->left, &errors, &bases);
      add_end(&r->right, &errors, &bases);
    }
    n += got;
  }
  fastq_batch_free(&batch);
  reads_reader_delete(reader);
  amplicon_analyzer_delete(a);
  if (n_reads != NULL)
    *n_reads = n;
  return bases == 0 ? 0.0 : (double)errors / (double)bases;
}

size_t auto_settings_decide_references(const int *segment,
    const size_t *votes, size_t n, bool *allowed) {
  size_t undecided = 0;
  for (size_t i = 0; i < n; i++)
    allowed[i] = true;
  for (int s = 1; s <= SEGMENTS_N; s++) {
    size_t n_refs = 0, total = 0, best = 0;
    for (size_t i = 0; i < n; i++) {
      if (segment[i] != s)
        continue;
      /* (of references with as many votes, the first) */
      if (n_refs == 0 || votes[i] > votes[best])
        best = i;
      n_refs++;
      total += votes[i];
    }
    if (n_refs < 2)
      continue;
    /* too few votes: no choice; without a clear lead (see
       auto_settings_clear_lead()) the reference with most votes is
       chosen all the same, as more reads would not change the shares */
    if (votes[best] < AUTO_REFERENCE_MIN_VOTES) {
      undecided++;
      continue;
    }
    for (size_t i = 0; i < n; i++)
      if (segment[i] == s)
        allowed[i] = i == best;
  }
  return undecided;
}

bool auto_settings_clear_lead(size_t best_votes, size_t segment_votes) {
  return (double)best_votes >=
         AUTO_REFERENCE_MIN_SHARE * (double)segment_votes;
}

/* add to <votes> those of the probes of the read <seq> (length <len>),
   see auto_settings_choose_references(); the probes which do not vote,
   as several references match them, are counted in <ambiguous>, by the
   segment <segment> of the reference found */
static void vote(reference_set_t *refs, const int *segment, const char *seq,
    size_t len, size_t *votes, size_t *ambiguous) {
  size_t probe_len = reference_set_probe_len(refs);
  if (len < probe_len)
    return;
  size_t room = len - probe_len, last = len;
  for (size_t j = 0; j < AUTO_REFERENCE_PROBES; j++) {
    size_t offset = room * (2 * j + 1) / (2 * AUTO_REFERENCE_PROBES);
    if (offset == last)
      continue;   /* (a short read: the same probe again) */
    last = offset;
    reference_hit_t hit;
    reference_set_locate_seeded(refs, seq + offset, probe_len, &hit);
    if (!hit.mapped)
      continue;
    if (hit.ambiguous)
      ambiguous[segment[hit.ref_index]]++;
    else
      votes[hit.ref_index]++;
  }
}

/* report the references chosen for <file> */
static void log_references(const trimmer_config_t *config, const char *file,
    const size_t *votes, const size_t *ambiguous, const bool *allowed,
    size_t n_reads, Logger *logger) {
  const fasta_t *fa = reference_index_fasta(config->index);
  for (int s = 1; s <= SEGMENTS_N; s++) {
    size_t n_refs = 0, n_allowed = 0, total = 0, chosen = 0;
    for (size_t i = 0; i < fa->n_records; i++) {
      if (reference_index_segment(config->index, i) != s)
        continue;
      n_refs++;
      total += votes[i];
      if (allowed[i]) {
        chosen = i;
        n_allowed++;
      }
    }
    if (n_refs < 2)
      continue;
    char name[16];
    segments_name(s, name, sizeof(name));
    if (n_allowed == 1)
      lprintf(logger, LOG_INFO, "reference of %s for %s: %s (%zu of %zu "
              "probes of the first %zu reads; %zu more match several "
              "references)%s\n", name, file, fa->records[chosen].name,
              votes[chosen], total, n_reads, ambiguous[s],
              auto_settings_clear_lead(votes[chosen], total) ? "" :
              "; about as close to other references");
    else
      lprintf(logger, LOG_INFO, "reference of %s for %s: not chosen, too "
              "few probes (%zu) in the first %zu reads (%zu more match "
              "several references): all %zu are used\n", name, file, total,
              n_reads, ambiguous[s], n_refs);
  }
}

bool* auto_settings_choose_references(const trimmer_config_t *config,
    const char *file, Logger *logger) {
  assert(config != NULL);
  assert(file != NULL);
  if (config->index == NULL)
    return NULL;
  const fasta_t *fa = reference_index_fasta(config->index);
  size_t n_refs = fa->n_records;
  /* a segment with several references is needed for a choice */
  size_t per_segment[SEGMENTS_N + 1] = {0};
  bool several = false;
  int *segment;
  XMALLOC(segment, sizeof(int) * (n_refs > 0 ? n_refs : 1));
  for (size_t i = 0; i < n_refs; i++) {
    segment[i] = reference_index_segment(config->index, i);
    if (segment[i] > 0 && ++per_segment[segment[i]] > 1)
      several = true;
  }
  if (!several) {
    free(segment);
    return NULL;
  }
  reference_set_t *refs = reference_set_new(config->index, config->probe_len,
      config->probe_error_rate);
  size_t *votes, ambiguous[SEGMENTS_N + 1] = {0};
  bool *allowed;
  XCALLOC(votes, n_refs, sizeof(size_t));
  XMALLOC(allowed, sizeof(bool) * n_refs);
  reads_reader_t *reader = reads_reader_new(file, logger);
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  size_t n = 0, undecided = 0;
  const size_t stages[] = {AUTO_REFERENCE_READS, AUTO_REFERENCE_MORE_READS};
  for (size_t stage = 0; stage < 2; stage++) {
    size_t got = 0;
    while (n < stages[stage] &&
           (got = reads_reader_read_batch(reader, &batch,
                                          stages[stage] - n)) > 0) {
      for (size_t i = 0; i < got; i++) {
        const fastq_record_t *rec = batch.records + i;
        vote(refs, segment, batch.data.s + rec->seq_from, rec->seq_len,
             votes, ambiguous);
      }
      n += got;
    }
    undecided = auto_settings_decide_references(segment, votes, n_refs,
                                                allowed);
    /* decided, or no more reads */
    if (undecided == 0 || n < stages[stage])
      break;
  }
  fastq_batch_free(&batch);
  reads_reader_delete(reader);
  reference_set_delete(refs);
  log_references(config, file, votes, ambiguous, allowed, n, logger);
  bool all = true;
  for (size_t i = 0; i < n_refs; i++)
    if (!allowed[i])
      all = false;
  free(votes);
  free(segment);
  if (all) {
    free(allowed);
    return NULL;
  }
  return allowed;
}

void auto_settings_decide(auto_settings_t *self, trimmer_config_t *config,
    const char *file, Logger *logger) {
  assert(self != NULL);
  assert(config != NULL);
  assert(file != NULL);
  if (!self->max_errors.pending)
    return;
  report_setting_t *s = &self->max_errors;
  s->measured = auto_settings_error_rate(config, file,
      AUTO_MAX_ERRORS_READS, &s->n_reads, logger);
  s->value = auto_settings_max_errors(s->measured);
  config->max_errors = s->value;
  free(self->file);
  self->file = xstrdup(file);
  s->file = self->file;
  s->pending = false;
  lprintf(logger, LOG_INFO, "Error rate of the primer alignments in the "
          "first %zu reads of %s: %s; errors in a primer: %g%% of its "
          "length (-e auto)\n", s->n_reads, file,
          fraction_percent_text(s->measured).s, 100 * s->value);
}

void auto_settings_free(auto_settings_t *self) {
  if (self == NULL)
    return;
  free(self->file);
  self->file = NULL;
  self->max_errors.file = NULL;
}

void auto_settings_test(void) {
  /* the choice of the reference of each segment: segment 4 has three
     references, segment 6 two, segment 1 one; a reference without a
     segment */
  {
    const int segment[] = {4, 4, 4, 6, 6, 1, 0};
    bool allowed[7];
    /* a clear winner in each segment */
    const size_t clear[] = {3, 90, 7, 40, 0, 0, 0};
    test_eq_size(auto_settings_decide_references(segment, clear, 7,
                                                 allowed), 0);
    test_assert(!allowed[0] && allowed[1] && !allowed[2]);
    test_assert(allowed[3] && !allowed[4]);
    test_assert(allowed[5] && allowed[6]);
    /* too few votes in segment 6: it keeps its references */
    const size_t few[] = {3, 90, 7, AUTO_REFERENCE_MIN_VOTES - 1, 0, 0, 0};
    test_eq_size(auto_settings_decide_references(segment, few, 7, allowed),
                 1);
    test_assert(!allowed[0] && allowed[1] && !allowed[2]);
    test_assert(allowed[3] && allowed[4]);
    /* no clear lead in segment 4 (a tie: the first is chosen), no vote
       at all in segment 6 */
    const size_t tie[] = {50, 50, 0, 0, 0, 0, 0};
    test_eq_size(auto_settings_decide_references(segment, tie, 7, allowed),
                 1);
    test_assert(allowed[0] && !allowed[1] && !allowed[2]);
    test_assert(allowed[3] && allowed[4]);
    test_assert(!auto_settings_clear_lead(50, 100));
    test_assert(auto_settings_clear_lead(60, 100));
    /* just the share needed */
    const size_t share[] = {60, 40, 0, 30, 20, 0, 0};
    test_eq_size(auto_settings_decide_references(segment, share, 7,
                                                 allowed), 0);
    test_assert(allowed[0] && !allowed[1] && !allowed[2]);
    test_assert(allowed[3] && !allowed[4]);
  }
  test_eq(auto_settings_max_errors(0.01), AUTO_MAX_ERRORS_LOW, double,
          "%g");
  test_eq(auto_settings_max_errors(AUTO_MAX_ERRORS_HIGH_RATE),
          AUTO_MAX_ERRORS_HIGH, double, "%g");
  test_eq(auto_settings_max_errors(0.1), AUTO_MAX_ERRORS_HIGH, double,
          "%g");
}
