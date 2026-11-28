/* (c) 2026 Giorgio Gonnella */

#include "json-sample.h"
#include "qc-flags.h"
#include "qc-summary.h"
#include "virus.h"
#include "segments.h"
#include "strutils.h"
#include "subtype-reads.h"

/* the evidence <e> in the HA or NA (<which>), as an object: the leading
   and the second subtype with their reads, among the reads with several
   voting probes */
static void write_evidence(FILE *f, size_t which,
    const subtype_reads_evidence_t *e) {
  char letter = subtype_reads_letter(which);
  fprintf(f, "{\"leading\": \"%c%d\", \"second\": ", letter, e->leading);
  if (e->second != 0)
    fprintf(f, "\"%c%d\"", letter, e->second);
  else
    fputs("null", f);
  fprintf(f, ", \"reads\": %zu, \"leading_reads\": %zu, "
          "\"second_reads\": %zu, \"second_expected\": %.2f, "
          "\"second_excess\": %.2f, \"second_share\": %.4f, "
          "\"above_threshold\": %s}", e->reads, e->leading_reads,
          e->second_reads, e->second_expected, e->second_excess,
          e->second_share, e->above_threshold ? "true" : "false");
}

/* the member "subtype_reads": the reads of each HA and NA subtype among
   the first reads of the sample (subtype-reads.h); null where there is
   nothing to say */
static void write_subtype_reads(json_writer_t *w, const counter_t *c) {
  static const char *const keys[] = {"ha", "na"};
  subtype_reads_evidence_t e[SUBTYPE_READS_N_SEGMENTS];
  bool any = false;
  for (size_t which = 0; which < SUBTYPE_READS_N_SEGMENTS; which++) {
    subtype_reads_evidence(counter_subtype_reads(c), which, e + which);
    any = any || e[which].known;
  }
  json_member(w, "subtype_reads");
  if (!any) {
    fputs("null", w->f);
    return;
  }
  json_open_block(w, '{');
  for (size_t which = 0; which < SUBTYPE_READS_N_SEGMENTS; which++) {
    json_member(w, keys[which]);
    if (e[which].known)
      write_evidence(w->f, which, e + which);
    else
      fputs("null", w->f);
  }
  json_close_block(w, '}', false);
}

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
  cleavage_result_t site;
  report_cleavage(r, c, &site);
  qc_flags(counter_reads(c), &qc, true, &site, counter_subtype_reads(c),
           &flags);
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
  write_subtype_reads(w, c);
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
