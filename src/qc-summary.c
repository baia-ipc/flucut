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

/* the depth of the reference <ref> and what it is made of */
static void measure(const fasta_t *refs, const counter_t *c, size_t ref,
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
  size_t *depth, *total;
  XMALLOC(depth, sizeof(size_t) * (len > 0 ? len : 1));
  XCALLOC(total, len > 0 ? len : 1, sizeof(size_t));
  double sum = 0, internal = 0, deletion = 0;
  for (size_t cat = 0; cat < COVER_N_CATEGORIES; cat++) {
    counter_coverage(c, ref, cat, depth);
    double cat_sum = 0;
    for (size_t p = 0; p < len; p++) {
      total[p] += depth[p];
      cat_sum += (double)depth[p];
    }
    sum += cat_sum;
    if (cat == COVER_INTERNAL_PLUS || cat == COVER_INTERNAL_MINUS)
      internal += cat_sum;
    else if (cat == COVER_DELETION)
      deletion += cat_sum;
  }
  size_t covered = 0, well = 0;
  for (size_t p = 0; p < len; p++) {
    if (total[p] > 0)
      covered++;
    if (total[p] >= QC_MIN_DEPTH)
      well++;
  }
  free(depth);
  free(total);
  if (len > 0) {
    q->mean_depth = sum / (double)len;
    q->covered = (double)covered / (double)len;
    q->well_covered = (double)well / (double)len;
  }
  if (sum > 0) {
    q->internal = internal / sum;
    q->deletion = deletion / sum;
  }
  report_lengths_t lengths;
  report_lengths(refs, c, ref, &lengths);
  q->complete = lengths.summary.total;
  q->short_complete = lengths.summary.short_count;
}

void qc_summary(const fasta_t *refs, const counter_t *c, qc_summary_t *out) {
  assert(c != NULL);
  assert(out != NULL);
  memset(out, 0, sizeof(qc_summary_t));
  if (refs == NULL)
    return;
  for (size_t r = 0; r < refs->n_records; r++) {
    const char *description = refs->records[r].description;
    int s = description == NULL ? 0 : segments_of_description(description);
    if (s == 0 || !counter_has_coverage(c, r))
      continue;
    qc_segment_t q;
    measure(refs, c, r, &q);
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
