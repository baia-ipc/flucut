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
  flag->severe = false;
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

/* the cleavage site of the haemagglutinin, when the reads do not settle
   it: the reason is what the reader needs, not the motif */
static void add_cleavage(qc_flags_t *out, const cleavage_result_t *site) {
  if (site == NULL || site->searched == 0)
    return;
  if (site->verdict == CLEAVAGE_TOO_FEW_READS)
    add(out, "cleavage_site", "only %zu reads over the HA cleavage site",
        site->reads);
  else if (site->disagrees_with_reference)
    add(out, "cleavage_site", "the reads over the HA cleavage site are "
        "of more than one haemagglutinin");
  else if (site->control_weak)
    add(out, "cleavage_site", "the reads disagree at the HA cleavage "
        "site and beside it (%.0f%% there)",
        100 * site->control_support);
  else if (site->verdict != CLEAVAGE_FIRM)
    add(out, "cleavage_site", "%.0f%% of the reads agree on the HA "
        "cleavage site", 100 * site->support);
}

/* the reads of a second HA or NA subtype, where their excess is above
   the threshold: what was observed, not what it means */
static void add_second_subtypes(qc_flags_t *out,
    const subtype_reads_t *subtypes) {
  static const char *const keys[] = {"second_ha_subtype",
                                     "second_na_subtype"};
  if (subtypes == NULL)
    return;
  for (size_t which = 0; which < SUBTYPE_READS_N_SEGMENTS; which++) {
    subtype_reads_evidence_t e;
    subtype_reads_evidence(subtypes, which, &e);
    if (e.above_threshold)
      add(out, keys[which], "reads of a second %s subtype (%c%d): %zu of "
          "%zu, %.1f%%", subtype_reads_gene(which),
          subtype_reads_letter(which), e.second, e.second_reads, e.reads,
          100.0 * (double)e.second_reads / (double)e.reads);
  }
}

double qc_flags_reads_placed(const counter_reads_t *s) {
  assert(s != NULL);
  size_t analyzed = counter_reads_analyzed(s);
  return analyzed == 0 ? 0 : (double)s->placed / (double)analyzed;
}

void qc_flags(const counter_reads_t *s, const qc_summary_t *qc, bool refs,
    const cleavage_result_t *site, const subtype_reads_t *subtypes,
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
  if (analyzed > 0 && qc_flags_reads_placed(s) < QC_FLAGS_PLACED) {
    add(out, "few_reads_placed", "only %.1f%% of the reads placed on a "
        "segment", 100.0 * qc_flags_reads_placed(s));
    out->flags[out->n - 1].severe = true;
  }
  add_segments(out, qc, "internal_priming", "of depth internally primed",
               internal_share, QC_FLAGS_INTERNAL);
  add_segments(out, qc, "internal_deletions",
               "of depth from reads with an internal deletion",
               deletion_share,
               QC_FLAGS_DELETION);
  add_cleavage(out, site);
  add_second_subtypes(out, subtypes);
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
  s.placed = 5000;
  s.origin[COUNTER_LEFT][ORIGIN_TERMINAL] = 5000;
  qc_flags(&s, &qc, true, NULL, NULL, &flags);
  test_eq_size(flags.n, 0);
  /* few reads, a quarter of them without a primer */
  s.reads_in = 400;
  s.n_ends[2] = 300;
  s.n_ends[0] = 100;
  qc_flags(&s, &qc, false, NULL, NULL, &flags);
  test_eq_size(flags.n, 2);
  test_eq_str(flags.flags[0].text, "few reads (400)");
  /* the cleavage site: flagged only when the reads do not settle it,
     and with the reason */
  cleavage_result_t site;
  memset(&site, 0, sizeof(site));
  s.reads_in = 5000;
  s.n_ends[2] = 5000;
  s.n_ends[0] = 0;
  site.searched = 900;
  site.reads = 600;
  site.support = 0.7;
  site.control_support = 0.9;
  site.verdict = CLEAVAGE_FIRM;
  qc_flags(&s, &qc, true, &site, NULL, &flags);
  test_eq_size(flags.n, 0);
  site.verdict = CLEAVAGE_CONFIRM;
  site.support = 0.3;
  qc_flags(&s, &qc, true, &site, NULL, &flags);
  test_eq_size(flags.n, 1);
  test_eq_str(flags.flags[0].key, "cleavage_site");
  test_eq_str(flags.flags[0].text,
              "30% of the reads agree on the HA cleavage site");
  site.control_weak = true;
  site.control_support = 0.42;
  qc_flags(&s, &qc, true, &site, NULL, &flags);
  test_eq_str(flags.flags[0].text, "the reads disagree at the HA cleavage "
              "site and beside it (42% there)");
  site.disagrees_with_reference = true;
  qc_flags(&s, &qc, true, &site, NULL, &flags);
  test_eq_str(flags.flags[0].text, "the reads over the HA cleavage site "
              "are of more than one haemagglutinin");
  site.verdict = CLEAVAGE_TOO_FEW_READS;
  site.reads = 7;
  qc_flags(&s, &qc, true, &site, NULL, &flags);
  test_eq_str(flags.flags[0].text, "only 7 reads over the HA cleavage site");
  /* nothing read, nothing said */
  memset(&site, 0, sizeof(site));
  qc_flags(&s, &qc, true, &site, NULL, &flags);
  test_eq_size(flags.n, 0);
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
  qc_flags(&s, &qc, true, NULL, NULL, &flags);
  test_eq_size(flags.n, 2);
  test_eq_str(flags.flags[0].text, "HA 25%, NA 12% of depth internally primed");
  test_eq_str(flags.flags[1].text, "PB2 90% of depth from reads "
              "with an internal deletion");
  /* the read ends not located */
  qc.not_located = 0.75;
  qc_flags(&s, &qc, true, NULL, NULL, &flags);
  test_eq_str(flags.flags[0].key, "not_located");
  /* reads of a second HA subtype above the threshold: 190 reads of H5
     (reference 0), 10 of H3 (reference 1), two probes each */
  const int segment[] = {4, 4}, number[] = {5, 3};
  const size_t h5[] = {0, 0}, h3[] = {1, 1};
  subtype_reads_t subtypes;
  subtype_reads_clear(&subtypes);
  subtypes.samples = 1;
  for (size_t i = 0; i < 200; i++)
    subtype_reads_add(&subtypes, segment, number, i < 190 ? h5 : h3, 2);
  qc.not_located = 0;
  memset(qc.segments, 0, sizeof(qc.segments));
  qc_flags(&s, &qc, true, NULL, &subtypes, &flags);
  test_eq_size(flags.n, 1);
  test_eq_str(flags.flags[0].key, "second_ha_subtype");
  test_eq_str(flags.flags[0].text,
              "reads of a second HA subtype (H3): 10 of 200, 5.0%");
  /* a single read of another subtype is not pointed out */
  subtype_reads_clear(&subtypes);
  subtypes.samples = 1;
  for (size_t i = 0; i < 200; i++)
    subtype_reads_add(&subtypes, segment, number, i < 199 ? h5 : h3, 2);
  qc_flags(&s, &qc, true, NULL, &subtypes, &flags);
  test_eq_size(flags.n, 0);
  /* few of the reads placed on a segment: the references describe a
     small part of the sample, a severe flag. At the threshold, or
     without references, nothing is said */
  memset(&qc, 0, sizeof(qc));
  s.placed = 1500;
  test_assert(qc_flags_reads_placed(&s) == 0.3);
  qc_flags(&s, &qc, true, NULL, NULL, &flags);
  test_eq_size(flags.n, 0);
  s.placed = 1499;
  qc_flags(&s, &qc, true, NULL, NULL, &flags);
  test_eq_size(flags.n, 1);
  s.placed = 70;
  qc_flags(&s, &qc, true, NULL, NULL, &flags);
  test_eq_size(flags.n, 1);
  test_eq_str(flags.flags[0].key, "few_reads_placed");
  test_eq_str(flags.flags[0].text,
              "only 1.4% of the reads placed on a segment");
  test_assert(flags.flags[0].severe);
  qc_flags(&s, &qc, false, NULL, NULL, &flags);
  test_eq_size(flags.n, 0);
  s.placed = 5000;
  qc.not_located = 0.75;
  qc_flags(&s, &qc, true, NULL, NULL, &flags);
  test_assert(!flags.flags[0].severe);
}
