/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include "qc-flags.h"
#include "test.h"

/* add a flag, with the text formatted as by printf */
static void add(qc_flags_t *out, const char *key, const char *fmt, ...)
  __attribute__((format(printf, 3, 4)));

static void add(qc_flags_t *out, const char *key, const char *fmt, ...) {
  if (out->n == QC_FLAGS_MAX)
    return;
  qc_flag_t *flag = out->flags + out->n++;
  flag->key = key;
  va_list ap;
  va_start(ap, fmt);
  vsnprintf(flag->text, sizeof(flag->text), fmt, ap);
  va_end(ap);
}

/* the share of the depth of a segment counted by a flag of the segments */
typedef double (*share_fn)(const qc_segment_t *q);

static double internal_share(const qc_segment_t *q) { return q->internal; }
static double deletion_share(const qc_segment_t *q) { return q->deletion; }

/* a flag listing the segments whose share is at least <min>, e.g.
   "HA 25%, NA 12% of depth internally primed" */
static void add_segments(qc_flags_t *out, const qc_summary_t *qc,
    const char *key, const char *what, share_fn share, double min) {
  char text[QC_FLAGS_TEXT_SIZE];
  size_t len = 0, n = 0;
  for (int s = 1; s <= SEGMENTS_N; s++) {
    const qc_segment_t *q = qc->segments + s;
    if (!qc_has_shares(q) || share(q) < min || len >= sizeof(text))
      continue;
    int written = snprintf(text + len, sizeof(text) - len, "%s%s %.0f%%",
                           n++ > 0 ? ", " : "", segments_gene(s),
                           100 * share(q));
    if (written > 0)
      len += (size_t)written;
  }
  if (n > 0)
    add(out, key, "%s %s", text, what);
}

void qc_flags(const counter_reads_t *s, const qc_summary_t *qc, bool refs,
    qc_flags_t *out) {
  assert(s != NULL);
  assert(qc != NULL);
  assert(out != NULL);
  out->n = 0;
  size_t analyzed = counter_reads_analyzed(s);
  if (s->reads_in < QC_FEW_READS)
    add(out, "few_reads", "few reads (%zu)", s->reads_in);
  if (analyzed > 0 &&
      (double)s->n_ends[0] >= QC_FLAGS_NO_PRIMER * (double)analyzed)
    add(out, "no_primer", "no primer in %.0f%% of the reads",
        100.0 * (double)s->n_ends[0] / (double)analyzed);
  if (!refs)
    return;
  if (qc->not_located >= QC_NOT_LOCATED)
    add(out, "not_located", "%.0f%% of the read ends not located in the "
        "references", 100 * qc->not_located);
  add_segments(out, qc, "internal_priming", "of depth internally primed",
               internal_share, QC_FLAGS_INTERNAL);
  add_segments(out, qc, "internal_deletions",
               "of depth from reads with an internal deletion",
               deletion_share,
               QC_FLAGS_DELETION);
}

void qc_flags_test(void) {
  counter_reads_t s;
  qc_summary_t qc;
  qc_flags_t flags;
  memset(&s, 0, sizeof(s));
  memset(&qc, 0, sizeof(qc));
  /* enough reads, all with primers and located: no flag */
  s.reads_in = 5000;
  s.n_ends[2] = 5000;
  s.origin[COUNTER_LEFT][ORIGIN_TERMINAL] = 5000;
  qc_flags(&s, &qc, true, &flags);
  test_eq_size(flags.n, 0);
  /* few reads, a quarter of them without a primer */
  s.reads_in = 400;
  s.n_ends[2] = 300;
  s.n_ends[0] = 100;
  qc_flags(&s, &qc, false, &flags);
  test_eq_size(flags.n, 2);
  test_eq_str(flags.flags[0].text, "few reads (400)");
  test_eq_str(flags.flags[1].text, "no primer in 25% of the reads");
  /* the segments above the thresholds are listed in one flag */
  s.reads_in = 5000;
  s.n_ends[2] = 5000;
  s.n_ends[0] = 0;
  qc.segments[4].found = qc.segments[6].found = true;
  qc.segments[1].found = true;
  qc.segments[1].mean_depth = qc.segments[4].mean_depth = 100;
  qc.segments[6].mean_depth = 100;
  qc.segments[4].internal = 0.248;
  qc.segments[6].internal = 0.12;
  qc.segments[1].internal = 0.05;
  qc.segments[1].deletion = 0.9;
  qc_flags(&s, &qc, true, &flags);
  test_eq_size(flags.n, 2);
  test_eq_str(flags.flags[0].text, "HA 25%, NA 12% of depth internally primed");
  test_eq_str(flags.flags[1].text, "PB2 90% of depth from reads "
              "with an internal deletion");
  /* the read ends not located */
  qc.not_located = 0.75;
  qc_flags(&s, &qc, true, &flags);
  test_eq_str(flags.flags[0].key, "not_located");
}
