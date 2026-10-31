/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdlib.h>
#include "alloc.h"
#include "amplicon.h"
#include "auto-settings.h"
#include "fastq-reader.h"
#include "reads-reader.h"
#include "reference.h"
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

int auto_settings_detect_subtype(const trimmer_config_t *config,
    const char *file, Logger *logger) {
  assert(config != NULL);
  assert(file != NULL);
  if (config->index == NULL)
    return -1;
  size_t ns = reference_index_n_subtypes(config->index);
  if (ns < 2)
    return -1;   /* one subtype or none: nothing to tell apart */
  amplicon_analyzer_t *a = amplicon_analyzer_new(config->scheme,
      config->max_errors, config->policy, config->window);
  reference_set_t *refs = reference_set_new(config->index, config->probe_len,
      config->probe_error_rate);
  reference_set_set_seeds_only(refs, config->seeds_only);
  reference_set_set_exhaustive(refs, config->exhaustive);
  amplicon_analyzer_use_references(a, refs, config->terminal_window);
  size_t *votes;
  XCALLOC(votes, ns, sizeof(size_t));
  reads_reader_t *reader = reads_reader_new(file, logger);
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  size_t n = 0, got, located = 0;
  while (n < AUTO_SUBTYPE_READS &&
         (got = reads_reader_read_batch(reader, &batch,
                                        AUTO_SUBTYPE_READS - n)) > 0) {
    for (size_t i = 0; i < got; i++) {
      const fastq_record_t *rec = batch.records + i;
      const amplicon_t *r = amplicon_analyzer_run(a,
          batch.data.s + rec->seq_from, rec->seq_len);
      const amplicon_end_t *e =
          r->left.found && r->left.ref.mapped ? &r->left :
          r->right.found && r->right.ref.mapped ? &r->right : NULL;
      if (e != NULL) {
        int st = reference_index_subtype(config->index, e->ref.ref_index);
        if (st >= 0) {
          votes[(size_t)st]++;
          located++;
        }
      }
    }
    n += got;
  }
  fastq_batch_free(&batch);
  reads_reader_delete(reader);
  amplicon_analyzer_delete(a);
  reference_set_delete(refs);
  int best = -1;
  size_t best_votes = 0;
  for (size_t j = 0; j < ns; j++)
    if (votes[j] > best_votes) {
      best_votes = votes[j];
      best = (int)j;
    }
  if (best >= 0)
    lprintf(logger, LOG_INFO, "subtype of %s: %s (%zu of %zu located "
            "reads)\n", file,
            reference_index_subtype_name(config->index, best), best_votes,
            located);
  free(votes);
  return best;
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
          "first %zu reads of %s: %.1f%%; errors in a primer: %g of its "
          "length (-e auto)\n", s->n_reads, file, 100 * s->measured,
          s->value);
}

void auto_settings_free(auto_settings_t *self) {
  if (self == NULL)
    return;
  free(self->file);
  self->file = NULL;
  self->max_errors.file = NULL;
}

void auto_settings_test(void) {
  test_eq(auto_settings_max_errors(0.01), AUTO_MAX_ERRORS_LOW, double,
          "%g");
  test_eq(auto_settings_max_errors(AUTO_MAX_ERRORS_HIGH_RATE),
          AUTO_MAX_ERRORS_HIGH, double, "%g");
  test_eq(auto_settings_max_errors(0.1), AUTO_MAX_ERRORS_HIGH, double,
          "%g");
}
