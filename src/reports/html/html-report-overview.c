/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "html-page.h"
#include "html-report-markup.h"
#include "html-report-overview.h"
#include "html-report-plots.h"
#include "html-report-sample.h"
#include "numutils.h"
#include "qc-flags.h"
#include "segments.h"

#define N_STATES (QC_SAMPLE_NOT_LOCATED + 1)

/* the samples in each state, as badges */
static void write_states(html_writer_t *w, const report_t *r) {
  static const markup_state_t badge[N_STATES] = {
    MARKUP_GOOD, MARKUP_WARNING, MARKUP_CRITICAL, MARKUP_NEUTRAL,
    MARKUP_NEUTRAL};
  static const char *const name[N_STATES] = {
    "complete", "partial", "with segments missing", "few reads",
    "not matched to the references"};
  size_t n[N_STATES] = {0}, waiting = 0;
  bool growing = r->running && r->options->live;
  for (size_t i = 0; i < r->n_sections; i++) {
    if (r->sections[i].counts == NULL) {
      waiting++;
      continue;
    }
    qc_summary_t qc;
    qc_summary(r->references, r->sections[i].counts, r->options->segment,
               &qc);
    n[qc.status]++;
  }
  markup_open(w, "p", "chips");
  markup_textf(w, "b", "count", "%zu samples", r->n_sections);
  for (size_t k = 0; k < N_STATES; k++)
    if (n[k] > 0)
      markup_badge(w, growing && k != QC_SAMPLE_COMPLETE ? MARKUP_NEUTRAL :
                   badge[k], "%zu %s", n[k], name[k]);
  if (waiting > 0)
    markup_badge(w, MARKUP_NEUTRAL, "%zu waiting", waiting);
  html_close(w);
}

/* the subtypes of the samples of the run, as a tile: the most frequent
   one, and how many samples have each (nothing if no sample has one) */
static void subtypes_tile(html_writer_t *w, const report_t *r) {
  enum { MAX_SUBTYPES = 6 };
  char names[MAX_SUBTYPES][QC_SUBTYPE_SIZE], list[160] = "";
  size_t counts[MAX_SUBTYPES] = {0}, n = 0, best = 0, len = 0;
  for (size_t i = 0; i < r->n_sections; i++) {
    if (r->sections[i].counts == NULL)
      continue;
    qc_summary_t qc;
    qc_summary(r->references, r->sections[i].counts, r->options->segment,
               &qc);
    if (qc.subtype[0] == '\0')
      continue;
    size_t k = 0;
    while (k < n && strcmp(names[k], qc.subtype) != 0)
      k++;
    if (k == MAX_SUBTYPES)
      continue;
    if (k == n)
      snprintf(names[n++], QC_SUBTYPE_SIZE, "%s", qc.subtype);
    counts[k]++;
  }
  if (n == 0)
    return;
  for (size_t k = 0; k < n; k++) {
    if (counts[k] > counts[best])
      best = k;
    int written = snprintf(list + len, sizeof(list) - len, "%s%s in %zu",
                           k > 0 ? ", " : "", names[k], counts[k]);
    if (written > 0 && len + (size_t)written < sizeof(list))
      len += (size_t)written;
  }
  markup_tile(w, n == 1 ? "subtype (closest references)" :
              "most frequent subtype", list, "%s", names[best]);
}

/* the resources of the run: while it is in progress, the threads busy
   now (since the previous update), the memory now and the elapsed time;
   at its end, the threads busy on average, the memory peak, the CPU
   time and the elapsed time */
static void resource_tiles(html_writer_t *w, const report_t *r) {
  const run_resources_t *u = r->resources;
  char sub[64];
  if (r->running) {
    snprintf(sub, sizeof(sub), "busy now: %.1f", u->load);
    markup_tile(w, "threads", sub, "%zu", u->threads);
    snprintf(sub, sizeof(sub), "peak %s",
             run_memory_text(u->memory_peak).s);
    markup_tile(w, "memory", sub, "%s",
                run_memory_text(u->memory > 0 ? u->memory :
                                u->memory_peak).s);
    snprintf(sub, sizeof(sub), "CPU time %s", run_duration_text(u->cpu).s);
    markup_tile(w, "elapsed", sub, "%s", run_duration_text(u->elapsed).s);
  } else {
    snprintf(sub, sizeof(sub), "busy on average: %.1f", u->load);
    markup_tile(w, "threads", sub, "%zu", u->threads);
    markup_tile(w, "memory peak", NULL, "%s",
                run_memory_text(u->memory_peak).s);
    markup_tile(w, "CPU time", NULL, "%s", run_duration_text(u->cpu).s);
    markup_tile(w, "elapsed", NULL, "%s", run_duration_text(u->elapsed).s);
  }
}

/* the throughput of the run */
static void write_tiles(html_writer_t *w, const report_t *r) {
  const counter_reads_t *s = counter_reads(r->total);
  size_t analyzed = counter_reads_analyzed(s), done = 0;
  char sub[128];
  for (size_t i = 0; i < r->n_sections; i++)
    if (r->sections[i].counts != NULL)
      done++;
  markup_open(w, "div", "tiles");
  if (r->running)
    markup_tile(w, "samples processed", NULL, "%zu/%zu", done,
                r->n_sections);
  else
    markup_tile(w, "samples", NULL, "%zu", r->n_sections);
  if (r->references != NULL)
    subtypes_tile(w, r);
  sample_reads_tile(w, s, r->resources);
  sample_length_tile(w, r->total);
  sample_written_tile(w, s);
  snprintf(sub, sizeof(sub), "one %s, none %s",
           percent_text(s->n_ends[1], analyzed).s,
           percent_text(s->n_ends[0], analyzed).s);
  markup_tile(w, "both primers", sub, "%s",
              percent_text(s->n_ends[2], analyzed).s);
  if (r->resources != NULL)
    resource_tiles(w, r);
  html_close(w);
}

/* the cell of a segment: the share of it which is well covered;
   <judged>: tinted by the state of the segment */
static void segment_cell(html_writer_t *w, const qc_summary_t *qc, int s,
    bool judged) {
  const qc_segment_t *q = qc->segments + s;
  html_open(w, "td");
  if (!qc_expects(qc, s)) {
    html_close(w);
    return;
  }
  if (!q->found) {
    html_attr(w, "class", judged ? "sc bad" : "sc");
    html_attr(w, "title", "%s: not found", segments_gene(s));
    html_raw(w, "&ndash;");
    html_close(w);
    return;
  }
  html_attr(w, "class", "sc %s", !judged ? "" :
            q->status == QC_SEGMENT_COMPLETE ? "ok" :
            q->status == QC_SEGMENT_PARTIAL ? "warn" : "bad");
  html_attr(w, "title", "%s: %s at \xe2\x89\xa5%d\xc3\x97, median depth "
            "%s\xc3\x97 (of full templates %s\xc3\x97)", segments_gene(s),
            fraction_percent_text(q->well_covered).s, QC_MIN_DEPTH,
            html_num(q->median_depth).s, html_num(q->median_full_depth).s);
  /* (rounded down: 100 only for the whole segment; a little more, as
     e.g. 100 * 0.29 is slightly less than 29 in floating point) */
  html_textf(w, "%d", (int)(100 * q->well_covered + 1e-9));
  html_close(w);
}

static double internal_share(const qc_segment_t *q) { return q->internal; }
static double deletion_share(const qc_segment_t *q) { return q->deletion; }

/* the cell of the segment with the greatest share of its depth from a
   kind of reads (e.g. "HA 25%"), emphasized from <min>; all the segments
   are listed when the cell is pointed at */
static void greatest_share_cell(html_writer_t *w, const qc_summary_t *qc,
    double (*share)(const qc_segment_t *), double min) {
  char all[160];
  size_t len = 0;
  int best = 0;
  for (int s = 1; s <= SEGMENTS_N; s++) {
    const qc_segment_t *q = qc->segments + s;
    if (!qc_has_shares(q))
      continue;
    if (best == 0 || share(q) > share(qc->segments + best))
      best = s;
    int n = snprintf(all + len, sizeof(all) - len, "%s%s %.1f%%",
                     len > 0 ? ", " : "", segments_gene(s),
                     100 * share(q));
    if (n > 0 && len + (size_t)n < sizeof(all))
      len += (size_t)n;
  }
  if (best == 0 || share(qc->segments + best) < 0.005) {
    markup_td_none(w);
    return;
  }
  html_open(w, "td");
  if (share(qc->segments + best) >= min)
    html_attr(w, "class", "hi");
  html_attr(w, "title", "%s", all);
  html_textf(w, "%s %.0f%%", segments_gene(best),
             100 * share(qc->segments + best));
  html_close(w);
}

/* the cell of the reads written by origin: a small stacked bar */
static void origin_cell(html_writer_t *w, const counter_reads_t *s) {
  static const char *const cls[READ_N_CATEGORIES] = {
    "c-full", "c-int", "c-del", "c-unp", "c-oth"};
  static const char *const name[READ_N_CATEGORIES] = {
    "full template", "internally primed", "internal deletion",
    "not placed", "other"};
  char all[250];
  size_t len = 0;
  for (size_t k = 0; k < READ_N_CATEGORIES; k++) {
    int n = snprintf(all + len, sizeof(all) - len, "%s%s %s",
                     k > 0 ? ", " : "", name[k],
                     percent_text(s->written[k], s->reads_out).s);
    if (n > 0 && len + (size_t)n < sizeof(all))
      len += (size_t)n;
  }
  html_open(w, "td");
  html_attr(w, "title", "%s", all);
  markup_open(w, "div", "sbar mini");
  for (size_t k = 0; k < READ_N_CATEGORIES; k++) {
    if (s->written[k] == 0)
      continue;
    html_open(w, "i");
    html_attr(w, "class", "%s", cls[k]);
    html_attr(w, "style", "width:%.1f%%",
              html_percent(s->written[k], s->reads_out));
    html_close(w);
  }
  html_close(w);
  html_close(w);
}

/* the row of the sample <i>; <max_reads>: the reads of the largest one */
static void sample_row(html_writer_t *w, const report_t *r, size_t i,
    size_t max_reads) {
  const report_section_t *g = r->sections + i;
  bool refs = r->references != NULL;
  markup_open_row(w);
  html_open(w, "td");
  html_open(w, "a");
  html_attr(w, "href", "#s%zu", i);
  html_text(w, g->name);
  html_close(w);
  html_close(w);
  if (g->counts == NULL) {
    html_open(w, "td");
    html_attr(w, "class", "wait");
    html_attr(w, "colspan", "%d", refs ? 17 : 4);
    html_text(w, "waiting");
    html_close(w);
    html_close(w);
    return;
  }
  const counter_reads_t *s = counter_reads(g->counts);
  size_t analyzed = counter_reads_analyzed(s);
  qc_summary_t qc;
  qc_summary(r->references, g->counts, r->options->segment, &qc);
  if (refs) {
    html_open(w, "td");
    sample_badge(w, r, &qc, true);
    html_close(w);
    if (qc.subtype[0] != '\0')
      markup_td(w, "%s", sample_subtype(&qc).s);
    else
      markup_td_none(w);
  }
  html_open(w, "td");
  html_attr(w, "class", "rd");
  html_attr(w, "style", "--f:%.1f%%", html_percent(s->reads_in, max_reads));
  html_text(w, html_num(s->reads_in).s);
  html_close(w);
  read_stats_summary_t lengths;
  read_stats_summarize(counter_read_stats(g->counts), COUNTER_STATS_INPUT,
                       &lengths);
  markup_td(w, "%s", html_num(lengths.n50).s);
  if (lengths.has_quality)
    markup_td(w, "Q%.0f", lengths.median_quality - 0.5);
  else
    markup_td_none(w);
  if (refs) {
    bool judged = qc.status <= QC_SAMPLE_INCOMPLETE &&
                  !(r->running && r->options->live);
    for (int k = 1; k <= SEGMENTS_N; k++)
      segment_cell(w, &qc, k, judged);
    /* (nothing is emphasized in a sample which cannot be judged) */
    greatest_share_cell(w, &qc, internal_share,
                        judged ? QC_FLAGS_INTERNAL : 2);
    greatest_share_cell(w, &qc, deletion_share,
                        judged ? QC_FLAGS_DELETION : 2);
  }
  if (refs)
    origin_cell(w, s);
  markup_td(w, "%s", percent_text(s->n_ends[2], analyzed).s);
  html_close(w);
}

/* a header cell; <title>: its tooltip (or NULL); <rows>, <cols>: the
   rows and columns it spans */
static void header_cell(html_writer_t *w, const char *text,
    const char *title, int rows, int cols) {
  html_open(w, "th");
  if (rows > 1)
    html_attr(w, "rowspan", "%d", rows);
  if (cols > 1)
    html_attr(w, "colspan", "%d", cols);
  if (title != NULL)
    html_attr(w, "title", "%s", title);
  html_raw(w, text);
  html_close(w);
}

/* the header of the table of samples, in two rows: the coverage of the
   segments under one heading, the other columns over both rows */
static void samples_header(html_writer_t *w) {
  char complete[160], coverage[96];
  snprintf(complete, sizeof(complete), "segments with at least %.0f%% of "
           "their length at a depth of %d or more, of the segments "
           "expected", 100 * QC_COMPLETE_FRACTION, QC_MIN_DEPTH);
  snprintf(coverage, sizeof(coverage), "coverage: %% of the segment at "
           "\xe2\x89\xa5%d\xc3\x97", QC_MIN_DEPTH);
  static const char *const after[] = {
    "int. primed", "int. deletions", "reads by origin", "both primers"};
  html_open(w, "thead");
  markup_open_row(w);
  header_cell(w, "sample", NULL, 2, 1);
  header_cell(w, "complete segments", complete, 2, 1);
  header_cell(w, "subtype", NULL, 2, 1);
  header_cell(w, "reads", NULL, 2, 1);
  header_cell(w, "N50", NULL, 2, 1);
  header_cell(w, "quality", NULL, 2, 1);
  header_cell(w, coverage, NULL, 1, SEGMENTS_N);
  for (size_t k = 0; k < sizeof(after) / sizeof(*after); k++)
    header_cell(w, after[k], NULL, 2, 1);
  html_close(w);
  markup_open_row(w);
  for (int s = 1; s <= SEGMENTS_N; s++)
    header_cell(w, segments_gene(s), NULL, 1, 1);
  html_close(w);
  html_close(w);
}

static void write_table(html_writer_t *w, const report_t *r) {
  size_t max_reads = 0;
  for (size_t i = 0; i < r->n_sections; i++)
    if (r->sections[i].counts != NULL &&
        counter_reads(r->sections[i].counts)->reads_in > max_reads)
      max_reads = counter_reads(r->sections[i].counts)->reads_in;
  markup_open(w, "div", "scroll");
  markup_open(w, "table", "samples");
  if (r->references != NULL)
    samples_header(w);
  else
    markup_header_row(w, "sample", "reads", "N50", "quality",
                      "both primers", NULL);
  for (size_t i = 0; i < r->n_sections; i++)
    sample_row(w, r, i, max_reads);
  html_close(w);
  html_close(w);
  if (r->references == NULL)
    return;
  markup_open(w, "div", "lg");
  html_textf(w, "PB2\xe2\x80\xa6NS: %% of the segment at \xe2\x89\xa5%d"
             "\xc3\x97 ", QC_MIN_DEPTH);
  markup_swatch(w, "t-ok", "complete");
  markup_swatch(w, "t-warn", "partial");
  markup_swatch(w, "t-bad", "low depth or not found");
  char na[128];
  snprintf(na, sizeof(na), "not judged (fewer than %s reads, read ends "
           "mostly not located, or still being sequenced)",
           html_num(QC_FEW_READS).s);
  markup_swatch(w, "t-na", na);
  html_text(w, "Reads by origin: ");
  markup_swatch(w, "c-full", "full template");
  markup_swatch(w, "c-int", "internally primed");
  markup_swatch(w, "c-del", "internal deletion");
  markup_swatch(w, "c-unp", "not placed");
  markup_swatch(w, "c-oth", "other");
  html_close(w);
}

/* the median depth of each segment in each sample, as a plot */
static void write_depths(html_writer_t *w, const report_t *r) {
  size_t *depths, n = 0;
  const char **names;
  bool *judged;
  XMALLOC(judged, sizeof(bool) * (r->n_sections + 1));
  XMALLOC(depths, sizeof(size_t) * (r->n_sections + 1) * SEGMENTS_N);
  XMALLOC(names, sizeof(char *) * (r->n_sections + 1));
  for (size_t i = 0; i < r->n_sections; i++) {
    if (r->sections[i].counts == NULL)
      continue;
    qc_summary_t qc;
    qc_summary(r->references, r->sections[i].counts, r->options->segment,
               &qc);
    for (int s = 1; s <= SEGMENTS_N; s++)
      depths[n * SEGMENTS_N + (size_t)s - 1] =
        qc.segments[s].found ? qc.segments[s].median_depth : 0;
    judged[n] = qc.status <= QC_SAMPLE_INCOMPLETE;
    names[n++] = r->sections[i].name;
  }
  if (n > 0) {
    markup_open(w, "section", "card depths");
    html_tag_text(w, "h3", "Median depth by segment");
    plots_segment_depths(w, depths, names, judged, n,
                         QC_MIN_DEPTH);
    markup_note(w, "One dot for each sample, in the order of the table "
                "(hollow: few reads, or not matched to the references).");
    html_close(w);
  }
  free(depths);
  free(names);
  free(judged);
}

void overview_write(html_writer_t *w, const report_t *r) {
  assert(r != NULL);
  if (r->references != NULL)
    write_states(w, r);
  write_tiles(w, r);
  write_table(w, r);
  if (r->references != NULL)
    write_depths(w, r);
}
