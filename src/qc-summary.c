/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "qc-summary.h"
#include "reference.h"
#include "report.h"
#include "test.h"

void qc_subtype(const char *ha, const char *na, char *buf) {
  buf[0] = '\0';
  /* the H of "H5N1", and its N */
  const char *h_end = strchr(ha, 'N'), *n_start = strchr(na, 'N');
  if (h_end == NULL && n_start == NULL)
    return;
  snprintf(buf, QC_SUBTYPE_SIZE, "%.*s%s",
           h_end == NULL ? 2 : (int)(h_end - ha), h_end == NULL ? "H?" : ha,
           n_start == NULL ? "N?" : n_start);
}

static int size_cmp(const void *a, const void *b) {
  size_t x = *(const size_t *)a, y = *(const size_t *)b;
  return x < y ? -1 : x > y;
}

void qc_reference(const fasta_t *refs, const counter_t *c, size_t ref,
    qc_segment_t *q) {
  const fasta_record_t *record = refs->records + ref;
  size_t len = record->seq_len;
  memset(q, 0, sizeof(qc_segment_t));
  q->found = true;
  q->ref = ref;
  q->ref_len = len;
  if (!reference_description_subtype(record->description, q->subtype,
                                     QC_SUBTYPE_SIZE))
    q->subtype[0] = '\0';
  size_t *depth, *total, *full;
  XMALLOC(depth, sizeof(size_t) * (len > 0 ? len : 1));
  XCALLOC(total, len > 0 ? len : 1, sizeof(size_t));
  XCALLOC(full, len > 0 ? len : 1, sizeof(size_t));
  double sum = 0, internal = 0, deletion = 0, minus = 0, stranded = 0;
  for (size_t cat = 0; cat < COVER_N_CATEGORIES; cat++) {
    counter_coverage(c, ref, cat, depth);
    double cat_sum = 0;
    bool is_full = cat == COVER_TERMINAL_PLUS || cat == COVER_TERMINAL_MINUS;
    for (size_t p = 0; p < len; p++) {
      total[p] += depth[p];
      if (is_full)
        full[p] += depth[p];
      cat_sum += (double)depth[p];
    }
    sum += cat_sum;
    if (cat == COVER_INTERNAL_PLUS || cat == COVER_INTERNAL_MINUS)
      internal += cat_sum;
    else if (cat == COVER_DELETION)
      deletion += cat_sum;
    if (cat == COVER_TERMINAL_MINUS || cat == COVER_INTERNAL_MINUS)
      minus += cat_sum;
    if (cat <= COVER_INTERNAL_MINUS)
      stranded += cat_sum;
  }
  size_t covered = 0, well = 0;
  for (size_t p = 0; p < len; p++) {
    if (total[p] > 0)
      covered++;
    if (total[p] >= QC_MIN_DEPTH)
      well++;
  }
  if (len > 0) {
    qsort(total, len, sizeof(size_t), size_cmp);
    q->median_depth = total[len / 2];
    qsort(full, len, sizeof(size_t), size_cmp);
    q->median_full_depth = full[len / 2];
  }
  q->bases = sum;
  q->minus = stranded > 0 ? minus / stranded : 0;
  free(depth);
  free(total);
  free(full);
  if (len > 0) {
    q->mean_depth = sum / (double)len;
    q->covered = (double)covered / (double)len;
    q->well_covered = (double)well / (double)len;
  }
  q->status = q->well_covered >= QC_COMPLETE_FRACTION ?
              QC_SEGMENT_COMPLETE : q->well_covered < QC_LOW_FRACTION ?
              QC_SEGMENT_LOW : QC_SEGMENT_PARTIAL;
  if (sum > 0) {
    q->internal = internal / sum;
    q->deletion = deletion / sum;
  }
  report_lengths_t lengths;
  report_lengths(refs, c, ref, &lengths);
  q->complete = lengths.summary.total;
  q->short_complete = lengths.summary.short_count;
}

bool qc_has_shares(const qc_segment_t *q) {
  return q->found && q->mean_depth >= QC_MIN_SHARE_DEPTH;
}

bool qc_expects(const qc_summary_t *self, int segment) {
  return self->only_segment == 0 || self->only_segment == segment;
}

/* the state of the sample with the counts <reads>, from that of its
   segments expected */
static void set_status(qc_summary_t *out, const counter_reads_t *reads) {
  size_t n_missing = 0;
  for (int s = 1; s <= SEGMENTS_N; s++) {
    if (!qc_expects(out, s))
      continue;
    out->n_expected++;
    if (out->segments[s].status == QC_SEGMENT_COMPLETE)
      out->n_complete++;
    else if (!out->segments[s].found)
      n_missing++;
  }
  size_t ends = 0, unknown = 0;
  for (size_t e = 0; e < COUNTER_N_ENDS; e++) {
    unknown += reads->origin[e][ORIGIN_UNKNOWN];
    ends += reads->origin[e][ORIGIN_UNKNOWN] +
            reads->origin[e][ORIGIN_TERMINAL] +
            reads->origin[e][ORIGIN_INTERNAL];
  }
  out->not_located = ends > 0 ? (double)unknown / (double)ends : 0;
  out->status = out->n_complete == out->n_expected ? QC_SAMPLE_COMPLETE :
                reads->reads_in < QC_FEW_READS ? QC_SAMPLE_FEW_READS :
                out->not_located >= QC_NOT_LOCATED ? QC_SAMPLE_NOT_LOCATED :
                n_missing > 0 ? QC_SAMPLE_INCOMPLETE : QC_SAMPLE_PARTIAL;
}

void qc_summary(const fasta_t *refs, const counter_t *c, int only_segment,
    qc_summary_t *out) {
  assert(c != NULL);
  assert(out != NULL);
  memset(out, 0, sizeof(qc_summary_t));
  out->only_segment = only_segment;
  if (refs == NULL) {
    set_status(out, counter_reads(c));
    return;
  }
  for (size_t r = 0; r < refs->n_records; r++) {
    const char *description = refs->records[r].description;
    int s = description == NULL ? 0 : segments_of_description(description);
    if (s == 0 || !counter_has_coverage(c, r))
      continue;
    qc_segment_t q;
    qc_reference(refs, c, r, &q);
    /* (of several references of a segment, the one with most reads) */
    if (!out->segments[s].found ||
        q.mean_depth > out->segments[s].mean_depth)
      out->segments[s] = q;
  }
  const char *first = NULL;
  for (int s = 1; s <= SEGMENTS_N; s++) {
    const qc_segment_t *q = out->segments + s;
    if (!q->found)
      continue;
    out->n_found++;
    if (q->subtype[0] == '\0')
      continue;
    if (first == NULL)
      first = q->subtype;
    else if (strcmp(first, q->subtype) != 0)
      out->mixed = true;
  }
  /* (segment 4: HA; segment 6: NA) */
  qc_subtype(out->segments[4].found ? out->segments[4].subtype : "",
             out->segments[6].found ? out->segments[6].subtype : "",
             out->subtype);
  set_status(out, counter_reads(c));
}

void qc_summary_test(void) {
  char subtype[QC_SUBTYPE_SIZE];
  qc_subtype("H5N1", "H5N1", subtype);
  test_eq_str(subtype, "H5N1");
  /* the HA and the NA closest to references of different subtypes */
  qc_subtype("H5N1", "H7N9", subtype);
  test_eq_str(subtype, "H5N9");
  qc_subtype("H10N8", "H3N2", subtype);
  test_eq_str(subtype, "H10N2");
  /* one of the two not found, or its reference without a subtype */
  qc_subtype("H5N1", "", subtype);
  test_eq_str(subtype, "H5N?");
  qc_subtype("", "H3N2", subtype);
  test_eq_str(subtype, "H?N2");
  qc_subtype("", "", subtype);
  test_eq_str(subtype, "");
}
