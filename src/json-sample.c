/* (c) 2026 Giorgio Gonnella */

#include "json-sample.h"
#include "qc-flags.h"
#include "qc-summary.h"
#include "virus.h"
#include "segments.h"
#include "strutils.h"

void json_sample_write(json_writer_t *w, const report_t *r,
    const counter_t *c) {
  static const char *const states[] = {
    "complete", "partial", "incomplete", "few_reads", "not_located"};
  static const char *const segment_states[] = {
    "missing", "low", "partial", "complete"};
  FILE *f = w->f;
  json_member(w, "sample");
  if (r->references == NULL) {
    fputs("null", f);
    return;
  }
  qc_summary_t qc;
  qc_flags_t flags;
  qc_summary(r->references, c, r->options->segment, &qc);
  qc_flags(counter_reads(c), &qc, true, &flags);
  json_open_block(w, '{');
  json_member(w, "state");
  fprintf(f, "\"%s\"", states[qc.status]);
  json_member(w, "virus");
  fprintf(f, "\"%s\"", virus_type_name(qc.virus));
  json_member(w, "subtype");
  if (qc.subtype[0] != '\0')
    fputs_json_string(qc.subtype, f);
  else
    fputs("null", f);
  json_count_member(w, "segments_expected", qc.n_expected);
  json_count_member(w, "segments_found", qc.n_found);
  json_count_member(w, "segments_complete", qc.n_complete);
  json_share_member(w, "ends_not_located", qc.not_located);
  json_member(w, "flags");
  fputc('[', f);
  bool first = true;
  for (size_t i = 0; i < flags.n; i++) {
    json_sep(f, &first);
    fprintf(f, "\"%s\"", flags.flags[i].key);
  }
  fputc(']', f);
  json_member(w, "segments");
  json_open_block(w, '[');
  for (int s = 1; s <= SEGMENTS_N; s++) {
    const qc_segment_t *q = qc.segments + s;
    if (!qc_expects(&qc, s))
      continue;
    json_member(w, NULL);
    /* (the segment number of the type of the sample: influenza B
       numbers PB1 1 and PB2 2, see segments.h) */
    fprintf(f, "{\"segment\": %d, \"gene\": \"%s\", \"state\": \"%s\"",
            segments_number(s, qc.virus), segments_gene(s),
            segment_states[q->status]);
    if (q->found) {
      fputs(", \"reference\": ", f);
      fputs_json_string(r->references->records[q->ref].name, f);
      fprintf(f, ", \"length\": %zu, \"mean_depth\": %.1f, "
              "\"median_depth\": %zu, \"median_full_depth\": %zu, "
              "\"covered\": %.4f, "
              "\"well_covered\": %.4f, \"internal\": %.4f, "
              "\"deletion\": %.4f, \"minus_strand\": %.4f, "
              "\"complete_amplicons\": %zu, \"short_complete\": %zu",
              q->ref_len, q->mean_depth, q->median_depth,
              q->median_full_depth, q->covered,
              q->well_covered, q->internal, q->deletion, q->minus,
              q->complete, q->short_complete);
    }
    fputc('}', f);
  }
  json_close_block(w, ']', true);
  json_close_block(w, '}', false);
}
