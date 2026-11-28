/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdbool.h>
#include "html-report-markers.h"
#include "html-report-markup.h"
#include "marker-counts.h"
#include "numutils.h"

/* the state of the badge of a verdict */
static markup_state_t verdict_state(size_t verdict) {
  switch (verdict) {
    case MARKER_FIRM: return MARKUP_GOOD;
    case MARKER_CONFIRM: return MARKUP_WARNING;
    default: return MARKUP_CRITICAL;
  }
}

/* the value of a marker, or a dash where the reads support none */
static const char* value_text(const marker_result_t *m) {
  bool held = m->found && (m->verdict == MARKER_FIRM ||
                           m->verdict == MARKER_CONFIRM);
  return held ? m->value : "\xe2\x80\x93";
}

/* one row: the marker, what the reads carry at it, and how far it can
   be trusted */
static void write_row(html_writer_t *w, size_t marker,
    const marker_result_t *m) {
  const marker_def_t *def = marker_def(marker);
  markup_open_row(w);
  markup_td(w, "%s", def->name);
  markup_td(w, "%s", value_text(m));
  if (m->reference_value[0] != '\0')
    markup_td(w, "%s", m->reference_value);
  else
    markup_td_none(w);
  markup_td(w, "%s of %s", percent_text(m->value_reads, m->reads).s,
            html_num(m->reads).s);
  markup_open(w, "td", NULL);
  markup_badge(w, verdict_state(m->verdict), "%s",
               marker_verdict_name(m->verdict));
  if (m->attention)
    markup_badge(w, MARKUP_WARNING, "worth a look");
  html_close(w);
  html_close(w);   /* row */
}

/*
   The values of a marker which enough reads carry, where there is more
   than one: the sample holds both, which is what reading the position
   off a consensus would not show
*/
static void write_mixtures(html_writer_t *w, const marker_counts_t *counts) {
  bool any = false;
  for (size_t i = 0; i < MARKER_N_MARKERS; i++) {
    marker_result_t m;
    marker_result_of(counts, i, &m);
    if (m.searched == 0 || m.supported < 2)
      continue;
    if (!any) {
      markup_details_start(w, "the markers at which the reads carry more "
                           "than one value");
      html_open(w, "table");
      markup_header_row(w, "marker", "value", "reads", "% of reads", NULL);
      any = true;
    }
    for (size_t k = 0; k < m.n_top; k++) {
      if (m.top[k].share < MARKER_MIN_SHARE ||
          m.top[k].reads < MARKER_MIN_VALUE_READS)
        continue;
      markup_open_row(w);
      markup_td(w, "%s", marker_def(i)->name);
      markup_td(w, "%s%s", m.top[k].value, m.top[k].attention ? " *" : "");
      markup_td(w, "%s", html_num(m.top[k].reads).s);
      markup_td(w, "%s", percent_text(m.top[k].reads, m.reads).s);
      html_close(w);
    }
  }
  if (any) {
    html_close(w);
    markup_details_end(w);
  }
}

void markers_write(html_writer_t *w, const counter_t *c) {
  assert(w != NULL);
  assert(c != NULL);
  const marker_counts_t *counts = counter_markers(c);
  size_t n = 0;
  for (size_t i = 0; i < MARKER_N_MARKERS; i++)
    if (counts->markers[i].searched > 0)
      n++;
  if (n == 0)
    return;
  markup_open(w, "section", "card");
  html_tag_text(w, "h3", "Markers");
  html_open(w, "table");
  markup_header_row(w, "marker", "reads", "reference", "agreeing",
                    "", NULL);
  for (size_t i = 0; i < MARKER_N_MARKERS; i++) {
    marker_result_t m;
    marker_result_of(counts, i, &m);
    if (m.searched == 0)
      continue;
    write_row(w, i, &m);
  }
  html_close(w);   /* table */
  write_mixtures(w, counts);
  markup_details_start(w, "what each marker is");
  html_open(w, "table");
  markup_header_row(w, "marker", "protein", "positions", "numbering",
                    "what it means", NULL);
  for (size_t i = 0; i < MARKER_N_MARKERS; i++) {
    if (counts->markers[i].searched == 0)
      continue;
    const marker_def_t *def = marker_def(i);
    markup_open_row(w);
    markup_td(w, "%s", def->name);
    markup_td(w, "%s", def->protein);
    markup_td(w, "%s", def->positions);
    markup_td(w, "%s", def->numbering);
    markup_td(w, "%s", def->meaning);
    html_close(w);
  }
  html_close(w);
  markup_details_end(w);
  markup_note(w, "The residues are read from the reads beside a conserved "
              "peptide, not from a consensus: the share is the reads "
              "carrying that value, so a residue present in a part of the "
              "sample is visible here.");
  html_close(w);   /* card */
}

bool markers_on_segment(int segment) {
  for (size_t i = 0; i < MARKER_N_MARKERS; i++)
    if (segment >= 0 && marker_def(i)->segment == (size_t)segment)
      return true;
  return false;
}

void markers_cell(html_writer_t *w, const counter_t *c, int segment) {
  assert(w != NULL);
  assert(c != NULL);
  const marker_counts_t *counts = counter_markers(c);
  html_open(w, "td");
  bool any = false;
  for (size_t i = 0; i < MARKER_N_MARKERS; i++) {
    const marker_def_t *def = marker_def(i);
    if (segment < 0 || def->segment != (size_t)segment)
      continue;
    marker_result_t m;
    marker_result_of(counts, i, &m);
    if (m.searched == 0)
      continue;
    if (any)
      html_text(w, " ");
    html_open(w, "span");
    html_attr(w, "class", "motif-cell");
    html_textf(w, "%s", value_text(&m));
    html_close(w);
    if (m.attention)
      markup_badge(w, MARKUP_WARNING, "%s", def->positions);
    any = true;
  }
  if (!any)
    html_text(w, "\xe2\x80\x93");
  html_close(w);
}
