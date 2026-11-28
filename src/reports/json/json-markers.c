/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include "json-markers.h"
#include "strutils.h"
#include "marker-counts.h"

/* a string member, or null where the string is empty */
static void string_or_null(json_writer_t *w, const char *key,
    const char *value) {
  json_member(w, key);
  if (value[0] == '\0')
    fputs("null", w->f);
  else
    fputs_json_string(value, w->f);
}

static void write_marker(json_writer_t *w, size_t marker,
    const marker_result_t *m) {
  const marker_def_t *def = marker_def(marker);
  FILE *f = w->f;
  json_open_block(w, '{');
  json_member(w, "name");
  fputs_json_string(def->name, f);
  json_member(w, "protein");
  fputs_json_string(def->protein, f);
  json_count_member(w, "segment", def->segment);
  json_member(w, "positions");
  fputs_json_string(def->positions, f);
  json_member(w, "numbering");
  fputs_json_string(def->numbering, f);
  json_member(w, "meaning");
  fputs_json_string(def->meaning, f);
  json_count_member(w, "reads_searched", m->searched);
  json_count_member(w, "reads", m->reads);
  json_count_member(w, "distinct_values", m->distinct);
  json_count_member(w, "supported_values", m->supported);
  string_or_null(w, "value", m->found ? m->value : "");
  json_count_member(w, "value_reads", m->value_reads);
  json_share_member(w, "support", m->support);
  string_or_null(w, "reference_value", m->reference_value);
  json_member(w, "differs_from_reference");
  fputs(m->differs_from_reference ? "true" : "false", f);
  /* what the same reads carry where the anchor matched: how much they
     agree where nothing is difficult */
  json_count_member(w, "anchor_reads", m->anchor_reads);
  json_share_member(w, "anchor_support", m->anchor_support);
  json_member(w, "reads_disagree_at_the_anchor");
  fputs(m->anchor_weak ? "true" : "false", f);
  /* the anchor as the reads carry it and as the reference has it: an
     anchor which matched somewhere else shows here */
  string_or_null(w, "anchor_residues", m->anchor_value);
  string_or_null(w, "reference_anchor_residues", m->reference_anchor);
  json_member(w, "anchor_differs_from_reference");
  fputs(m->anchor_differs ? "true" : "false", f);
  json_member(w, "attention");
  fputs(m->attention ? "true" : "false", f);
  json_share_member(w, "attention_share", m->attention_share);
  json_member(w, "verdict");
  fprintf(f, "\"%s\"", marker_verdict_name(m->verdict));
  /* the values the reads carry, most frequent first: a value carried by
     a part of the reads is a part of the sample, not a call */
  json_member(w, "values");
  json_open_block(w, '[');
  for (size_t i = 0; i < m->n_top; i++) {
    json_member(w, NULL);
    fputs("{\"value\": ", f);
    fputs_json_string(m->top[i].value, f);
    fprintf(f, ", \"reads\": %zu, \"share\": %.4f, \"attention\": %s}",
            m->top[i].reads, m->top[i].share,
            m->top[i].attention ? "true" : "false");
  }
  json_close_block(w, ']', true);
  json_close_block(w, '}', true);
}

void json_markers_write(json_writer_t *w, const counter_t *c) {
  assert(w != NULL);
  assert(c != NULL);
  const marker_counts_t *counts = counter_markers(c);
  json_member(w, "markers");
  json_open_block(w, '[');
  for (size_t i = 0; i < MARKER_N_MARKERS; i++) {
    marker_result_t m;
    marker_result_of(counts, i, &m);
    if (m.searched == 0)
      continue;
    json_member(w, NULL);
    write_marker(w, i, &m);
  }
  json_close_block(w, ']', true);
}
