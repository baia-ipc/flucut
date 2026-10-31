/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "file.h"
#include "html-page.h"
#include "html-report.h"
#include "html-template.h"
#include "html-writer.h"
#include "segments.h"
#include "strbuf.h"

#include "report-logo-data.inc"

/* the static parts of the report (html-report-*.html), embedded by the
   Makefile */
static const char *const run_tmpl =
#include "html-report-run-html.inc"
  ;
static const char *const options_tmpl =
#include "html-report-options-html.inc"
  ;
static const char *const style_tmpl =
#include "html-report-style-html.inc"
  ;
static const char *const steps_note_tmpl =
#include "html-report-steps-note-html.inc"
  ;
static const char *const primers_note_tmpl =
#include "html-report-primers-note-html.inc"
  ;
static const char *const sites_note_tmpl =
#include "html-report-sites-note-html.inc"
  ;
static const char *const no_sites_note_tmpl =
#include "html-report-no-sites-note-html.inc"
  ;
static const char *const lengths_note_tmpl =
#include "html-report-lengths-note-html.inc"
  ;

/* height of the logo (pixels; the images are twice as large) */
#define REPORT_LOGO_HEIGHT 75
/* the logo in the top bar of the dashboard is smaller */
#define DASHBOARD_LOGO_HEIGHT 40

/* --- markup helpers (as in Lumipore) --- */

/* start an element of the class <cls> */
static void open_class(html_writer_t *w, const char *tag, const char *cls) {
  html_open(w, tag);
  html_attr(w, "class", "%s", cls);
}

/* a paragraph of the class "note", with the text formatted as by printf */
static void note(html_writer_t *w, const char *fmt, ...)
  __attribute__((format(printf, 2, 3)));

static void note(html_writer_t *w, const char *fmt, ...) {
  open_class(w, "p", "note");
  va_list ap;
  va_start(ap, fmt);
  html_vtextf(w, fmt, ap);
  va_end(ap);
  html_close(w);
}

/* start a table row (on one line) */
static void open_row(html_writer_t *w) {
  html_open(w, "tr");
  html_inline(w);
}

/* a header row: the cells contain the given markup (constant strings),
   up to a NULL */
static void header_row(html_writer_t *w, ...) {
  open_row(w);
  va_list ap;
  va_start(ap, w);
  for (const char *h = va_arg(ap, const char *); h != NULL;
       h = va_arg(ap, const char *)) {
    html_open(w, "th");
    html_raw(w, h);
    html_close(w);
  }
  va_end(ap);
  html_close(w);
}

/* a cell with the text formatted as by printf */
#define td(W, ...) html_tag_textf((W), "td", __VA_ARGS__)

/* a cell without value */
static void td_none(html_writer_t *w) {
  html_open(w, "td");
  html_raw(w, "&ndash;");
  html_close(w);
}

/* a tile with the value formatted as by printf */
static void tile(html_writer_t *w, const char *label, const char *fmt, ...)
  __attribute__((format(printf, 3, 4)));

static void tile(html_writer_t *w, const char *label, const char *fmt,
    ...) {
  char value[64];
  va_list ap;
  va_start(ap, fmt);
  vsnprintf(value, sizeof(value), fmt, ap);
  va_end(ap);
  html_tile(w, value, label);
}

/* a count, with its percentage of <total> */
static void td_count(html_writer_t *w, size_t n, size_t total) {
  td(w, "%s", html_num(n).s);
  td(w, "%.1f%%", html_percent(n, total));
}

/* --- sections --- */

static void write_tiles(html_writer_t *w, const report_t *r,
    const counter_reads_t *s) {
  tile(w, "reads", "%s", html_num(s->reads_in).s);
  tile(w, "written", "%.1f%%", html_percent(s->reads_out, s->reads_in));
  tile(w, "with both primers", "%.1f%%",
       html_percent(s->n_ends[2], s->reads_in));
  if (r->references != NULL)
    tile(w, "internally primed ends", "%s",
         html_num(s->origin[COUNTER_LEFT][ORIGIN_INTERNAL] +
                  s->origin[COUNTER_RIGHT][ORIGIN_INTERNAL]).s);
  else
    tile(w, "primer dimers", "%s", html_num(s->dimers).s);
}

/* the reads at each step */
static void write_steps(html_writer_t *w, const counter_reads_t *s) {
  const char *labels[] = {"all reads", "both primers", "one primer",
                          "no primer", "primer dimers",
                          "too short (not written)", "written",
                          "minus strand", "adapter at the start",
                          "adapter at the end"};
  size_t values[] = {s->reads_in, s->n_ends[2], s->n_ends[1], s->n_ends[0],
                     s->dimers, s->too_short, s->reads_out, s->minus,
                     s->adapter_start, s->adapter_end};
  const char *classes[] = {"p", "s", "s", "p", "p", "p", "s", "p", "p",
                           "p"};
  html_hbars(w, sizeof(values) / sizeof(values[0]), labels, values,
             classes, s->reads_in, 0);
  html_template(w, steps_note_tmpl);
}

/* the primers found at each end, by primer of the scheme */
static void write_primers(html_writer_t *w, const report_t *r,
    const counter_t *c) {
  size_t total = counter_reads(c)->reads_in;
  html_open(w, "table");
  header_row(w, "end", "primer", "reads", "% of the reads", NULL);
  for (size_t e = 0; e < COUNTER_N_ENDS; e++) {
    const primer_t *p = e == COUNTER_LEFT ? r->scheme->left :
                        r->scheme->right;
    for (size_t v = 0; v < counter_n_variants(c, e); v++) {
      open_row(w);
      td(w, "%s", report_end_name(e));
      td(w, "%s", p[v].name);
      td_count(w, counter_primer(c, e, v), total);
      html_close(w);
    }
  }
  html_close(w);
  html_template(w, primers_note_tmpl);
}

/* the read ends by origin, and the internal priming sites */
static void write_sites(html_writer_t *w, const report_t *r,
    const counter_t *c) {
  if (r->references == NULL) {
    html_template(w, no_sites_note_tmpl);
    return;
  }
  const counter_reads_t *s = counter_reads(c);
  html_open(w, "table");
  header_row(w, "end", "terminal", "internal", "not located", NULL);
  for (size_t e = 0; e < COUNTER_N_ENDS; e++) {
    open_row(w);
    td(w, "%s", report_end_name(e));
    td(w, "%s", html_num(s->origin[e][ORIGIN_TERMINAL]).s);
    td(w, "%s", html_num(s->origin[e][ORIGIN_INTERNAL]).s);
    td(w, "%s", html_num(s->origin[e][ORIGIN_UNKNOWN]).s);
    html_close(w);
  }
  html_close(w);
  size_t n;
  report_site_t *sites = report_sites(r->references, c, &n);
  if (n == 0) {
    note(w, "No read end was classified as internally primed.");
    return;
  }
  size_t shown = n < HTML_REPORT_MAX_SITES ? n : HTML_REPORT_MAX_SITES;
  html_open(w, "table");
  header_row(w, "reference", "end", "position", "range", "reads",
             "% of the reads", NULL);
  for (size_t i = 0; i < shown; i++) {
    const report_site_t *x = sites + i;
    open_row(w);
    const fasta_record_t *ref = r->references->records + x->ref;
    char label[REPORT_LABEL_SIZE];
    segments_label(ref->name, ref->description, label, sizeof(label));
    td(w, "%s", label);
    td(w, "%s", report_end_name(x->end));
    td(w, "%zu", x->pos);
    if (x->from == x->to)
      td_none(w);
    else
      td(w, "%zu-%zu", x->from, x->to);
    td_count(w, x->reads, s->reads_in);
    html_close(w);
  }
  html_close(w);
  char gap[32], n_shown[32], max[32];
  snprintf(gap, sizeof(gap), "%d", REPORT_SITE_GAP);
  snprintf(n_shown, sizeof(n_shown), "%zu of %zu", shown, n);
  snprintf(max, sizeof(max), "%d", HTML_REPORT_MAX_SITES);
  html_template(w, sites_note_tmpl, "gap", gap, "shown", n_shown,
                "max", max);
  free(sites);
}

/* a bin of the lengths, as "<from>-<to>" */
static void td_bin(html_writer_t *w, size_t bin) {
  if (bin == LENGTHS_N_BINS - 1)
    td(w, "%d+", LENGTHS_MAX);
  else
    td(w, "%zu-%zu", bin * LENGTHS_BIN,
       bin * LENGTHS_BIN + LENGTHS_BIN - 1);
}

/* the complete amplicons of each group of the lengths, and their
   histograms */
static void write_lengths(html_writer_t *w, const report_t *r,
    const counter_t *c) {
  const lengths_t *l = counter_lengths(c);
  size_t n_groups = lengths_n_groups(l);
  html_open(w, "table");
  header_row(w, r->references != NULL ? "reference" : "group",
             "length", "complete amplicons", "most at", "short",
             "short, most at", NULL);
  for (size_t g = 0; g < n_groups; g++) {
    report_lengths_t x;
    report_lengths(r->references, c, g, &x);
    if (x.summary.total == 0)
      continue;
    open_row(w);
    td(w, "%s", x.label);
    if (x.ref_len > 0)
      td(w, "%s", html_num(x.ref_len).s);
    else
      td_none(w);
    td(w, "%s", html_num(x.summary.total).s);
    td_bin(w, x.summary.main_bin);
    if (x.ref_len == 0) {
      td_none(w);
      td_none(w);
    } else {
      td(w, "%s (%.1f%%)", html_num(x.summary.short_count).s,
         html_percent(x.summary.short_count, x.summary.total));
      if (x.summary.short_count > 0)
        td_bin(w, x.summary.short_bin);
      else
        td_none(w);
    }
    html_close(w);
  }
  html_close(w);
  char percent[32], bin[32];
  snprintf(percent, sizeof(percent), "%.0f", 100 * LENGTHS_SHORT_FRACTION);
  snprintf(bin, sizeof(bin), "%d", LENGTHS_BIN);
  html_template(w, lengths_note_tmpl, "percent", percent, "bin", bin);
  /* the histograms, in a grid */
  html_tag_text(w, "h3", "Lengths of the complete amplicons");
  open_class(w, "div", "plots");
  size_t hist[LENGTHS_N_BINS];
  for (size_t g = 0; g < n_groups; g++) {
    report_lengths_t x;
    report_lengths(r->references, c, g, &x);
    if (x.summary.total == 0)
      continue;
    for (size_t b = 0; b < LENGTHS_N_BINS; b++)
      hist[b] = lengths_count(l, g, b);
    open_class(w, "div", "plot");
    html_tag_text(w, "b", x.label);
    html_histogram(w, hist, LENGTHS_N_BINS, LENGTHS_BIN,
                   "length between the primers", "bases",
                   "No complete amplicons.");
    html_close(w);
  }
  html_close(w);
}

/* the colors and labels of the coverage categories (COVER_*): the
   internally primed reads stand out (amber, blue), the full-length ones
   are green, the reads whose 5' origin is unknown are grey */
static const char *const cover_color[COVER_N_CATEGORIES] = {
  "#0b9e73", "#7fc7a9", "#e8a013", "#2e7fe0", "#9e9e9e", "#cccccc"};
static const char *const cover_label[COVER_N_CATEGORIES] = {
  "terminal +", "terminal \xe2\x88\x92", "internal +", "internal \xe2\x88\x92",
  "no 5' primer", "5' probe not located"};

/* the order in which the categories are stacked, from the bottom: the
   flat categories first, so that the bands stay horizontal. A stacked
   band is horizontal only when everything below it is flat; the terminal
   reads and the reads of unknown origin cover the whole segment evenly
   (flat), while the internally primed reads are absent before their site
   and present after it (a step). So the grey (unknown) and the green
   (terminal) are at the bottom as flat bands, and the internally primed
   reads (amber, blue) are on top, where their step up at the site is the
   artefact and the bands below stay horizontal */
static const int cover_order[COVER_N_CATEGORIES] = {
  COVER_UNMAPPED, COVER_NO_PRIMER, COVER_TERMINAL_PLUS, COVER_TERMINAL_MINUS,
  COVER_INTERNAL_PLUS, COVER_INTERNAL_MINUS};

/* a stacked-area coverage plot of one reference, by read category:
   x the position on the reference, y the number of located reads covering
   it, stacked by category (see counter_coverage()); the largest internal
   priming site is marked */
static void coverage_plot(html_writer_t *w, const counter_t *c, size_t ref,
    size_t len, const char *label) {
  enum { ML = 46, MR = 10, MT = 10, PH = 150, PW = 760 };
  size_t inner = PW - ML - MR;
  size_t nc = len < inner ? len : inner;
  if (nc < 2)
    return;
  size_t *cov;
  XMALLOC(cov, sizeof(size_t) * COVER_N_CATEGORIES * len);
  size_t ymax = 0;
  for (size_t k = 0; k < COVER_N_CATEGORIES; k++)
    counter_coverage(c, ref, k, cov + k * len);
  for (size_t p = 0; p < len; p++) {
    size_t sum = 0;
    for (size_t k = 0; k < COVER_N_CATEGORIES; k++)
      sum += cov[k * len + p];
    if (sum > ymax)
      ymax = sum;
  }
  if (ymax == 0) {
    free(cov);
    note(w, "No located reads.");
    return;
  }
  double *base;
  XCALLOC(base, nc, sizeof(double));
  html_open(w, "svg");
  html_attr(w, "viewBox", "0 0 %d %d", PW, MT + PH + 20);
  html_attr(w, "width", "100%%");
  html_attr(w, "style", "max-width:%dpx;height:auto", PW);
  html_attr(w, "role", "img");
  /* the axes */
  html_elem(w, "line");
  html_attr(w, "x1", "%d", ML); html_attr(w, "y1", "%d", MT);
  html_attr(w, "x2", "%d", ML); html_attr(w, "y2", "%d", MT + PH);
  html_attr(w, "stroke", "currentColor"); html_attr(w, "opacity", "0.4");
  html_elem(w, "line");
  html_attr(w, "x1", "%d", ML); html_attr(w, "y1", "%d", MT + PH);
  html_attr(w, "x2", "%zu", ML + inner); html_attr(w, "y2", "%d", MT + PH);
  html_attr(w, "stroke", "currentColor"); html_attr(w, "opacity", "0.4");
  /* the stacked areas */
  strbuf_t pts;
  strbuf_init(&pts);
  for (size_t ki = 0; ki < COVER_N_CATEGORIES; ki++) {
    size_t k = (size_t)cover_order[ki];
    const size_t *v = cov + k * len;
    strbuf_clear(&pts);
    char buf[32];
    for (size_t i = 0; i < nc; i++) {
      size_t p = (size_t)((double)i * (len - 1) / (nc - 1));
      double x = ML + (double)i * inner / (nc - 1);
      double y = MT + PH - PH * (base[i] + v[p]) / ymax;
      int n = snprintf(buf, sizeof(buf), "%s%.1f,%.1f", i ? " " : "", x, y);
      strbuf_append(&pts, buf, (size_t)n);
    }
    for (size_t ii = nc; ii-- > 0;) {
      double x = ML + (double)ii * inner / (nc - 1);
      double y = MT + PH - PH * base[ii] / ymax;
      int n = snprintf(buf, sizeof(buf), " %.1f,%.1f", x, y);
      strbuf_append(&pts, buf, (size_t)n);
    }
    html_elem(w, "polygon");
    html_attr(w, "points", "%.*s", (int)pts.len, pts.s);
    html_attr(w, "fill", "%s", cover_color[k]);
    for (size_t i = 0; i < nc; i++) {
      size_t p = (size_t)((double)i * (len - 1) / (nc - 1));
      base[i] += v[p];
    }
  }
  strbuf_free(&pts);
  /* the largest internal priming site of the reference */
  const size_t *il = counter_sites(c, ref, COUNTER_LEFT);
  const size_t *ir = counter_sites(c, ref, COUNTER_RIGHT);
  size_t best = 0, bestp = 0;
  for (size_t p = 0; p < len; p++) {
    size_t v = (il ? il[p] : 0) + (ir ? ir[p] : 0);
    if (v > best) { best = v; bestp = p; }
  }
  if (best > 0) {
    double x = ML + (double)bestp * inner / (len > 1 ? len - 1 : 1);
    html_elem(w, "line");
    html_attr(w, "x1", "%.1f", x); html_attr(w, "y1", "%d", MT);
    html_attr(w, "x2", "%.1f", x); html_attr(w, "y2", "%d", MT + PH);
    html_attr(w, "stroke", "#d7263d"); html_attr(w, "stroke-width", "1");
  }
  /* the labels */
  html_open(w, "text");
  html_attr(w, "x", "%d", ML - 4); html_attr(w, "y", "%d", MT + 4);
  html_attr(w, "text-anchor", "end"); html_attr(w, "font-size", "10");
  html_attr(w, "fill", "currentColor");
  html_inline(w); html_textf(w, "%zu", ymax); html_close(w);
  html_open(w, "text");
  html_attr(w, "x", "%zu", ML + inner); html_attr(w, "y", "%d", MT + PH + 14);
  html_attr(w, "text-anchor", "end"); html_attr(w, "font-size", "10");
  html_attr(w, "fill", "currentColor");
  html_inline(w); html_textf(w, "%zu", len); html_close(w);
  html_close(w);   /* svg */
  (void)label;
  free(base);
  free(cov);
}

/* the coverage plots of the references, in a grid, with a legend (see
   coverage_plot()); only the references with located reads */
static void write_coverage(html_writer_t *w, const report_t *r,
    const counter_t *c) {
  if (r->references == NULL)
    return;
  size_t shown = 0;
  for (size_t g = 0; g < r->references->n_records; g++)
    if (counter_has_coverage(c, g))
      shown++;
  if (shown == 0) {
    note(w, "No reads located in the references.");
    return;
  }
  open_class(w, "div", "coverlegend");
  for (size_t k = 0; k < COVER_N_CATEGORIES; k++) {
    html_open(w, "span");
    html_open(w, "i");
    html_attr(w, "style", "background:%s", cover_color[k]);
    html_close(w);
    html_text(w, cover_label[k]);
    html_close(w);
  }
  html_close(w);
  /* the references in order of segment number (1-8), those without a
     recognized number last, keeping the file order among ties */
  size_t n = r->references->n_records;
  size_t *order;
  XMALLOC(order, sizeof(size_t) * n);
  for (size_t i = 0; i < n; i++)
    order[i] = i;
  for (size_t i = 1; i < n; i++) {   /* a stable insertion sort */
    size_t gi = order[i];
    int si = segments_of_description(r->references->records[gi].description);
    if (si == 0) si = 999;
    size_t j = i;
    while (j > 0) {
      int sj = segments_of_description(
          r->references->records[order[j - 1]].description);
      if (sj == 0) sj = 999;
      if (sj <= si)
        break;
      order[j] = order[j - 1];
      j--;
    }
    order[j] = gi;
  }
  for (size_t i = 0; i < n; i++) {
    size_t g = order[i];
    if (!counter_has_coverage(c, g))
      continue;
    const fasta_record_t *ref = r->references->records + g;
    char label[REPORT_LABEL_SIZE];
    segments_label(ref->name, ref->description, label, sizeof(label));
    open_class(w, "div", "coverplot");
    html_tag_text(w, "b", label);
    coverage_plot(w, c, g, ref->seq_len, label);
    html_close(w);
  }
  free(order);
}

/* a compact line of basic run statistics (platform, bases and mean read
   lengths, minus strand), below the tiles; the data are collected anyway */
static void write_stats(html_writer_t *w, const report_t *r,
    const counter_reads_t *s) {
  open_class(w, "p", "stats");
  html_textf(w, "Platform: %s. ",
             r->options->illumina ? "Illumina (short reads)" :
             "Nanopore (long reads)");
  html_textf(w, "Input: %s reads, %s bases",
             html_num(s->reads_in).s, html_num(s->bases_in).s);
  if (s->reads_in > 0)
    html_textf(w, ", mean %s bp",
               html_num(s->bases_in / s->reads_in).s);
  html_textf(w, ". Written: %s reads, %s bases",
             html_num(s->reads_out).s, html_num(s->bases_out).s);
  if (s->reads_out > 0)
    html_textf(w, ", mean %s bp",
               html_num(s->bases_out / s->reads_out).s);
  html_textf(w, ". Minus strand: %s (%.1f%%).",
             html_num(s->minus).s, html_percent(s->minus, s->reads_in));
  html_close(w);
}

static void write_section(html_writer_t *w, const report_t *r,
    const counter_t *c) {
  const counter_reads_t *s = counter_reads(c);
  open_class(w, "div", "tiles");
  write_tiles(w, r, s);
  html_close(w);
  write_stats(w, r, s);
  if (s->reads_in == 0) {
    note(w, "No reads.");
    return;
  }
  html_tag_text(w, "h3", "Reads at each step");
  write_steps(w, s);
  html_tag_text(w, "h3", "Primers found");
  write_primers(w, r, c);
  html_tag_text(w, "h3", "Internal priming");
  write_sites(w, r, c);
  html_tag_text(w, "h3", "Complete amplicons");
  if (s->n_ends[2] == 0)
    note(w, "No complete amplicons.");
  else
    write_lengths(w, r, c);
  if (r->references != NULL) {
    html_tag_text(w, "h3", "Coverage by read origin");
    write_coverage(w, r, c);
  }
}

/* the groups: reads, written, with both primers (and the internal ends);
   the group names link to the sections of the report */
static void write_groups_table(html_writer_t *w, const report_t *r) {
  html_open(w, "table");
  if (r->references != NULL)
    header_row(w, "group", "files", "reads", "written", "both primers",
               "internal ends", NULL);
  else
    header_row(w, "group", "files", "reads", "written", "both primers",
               NULL);
  for (size_t i = 0; i < r->n_sections; i++) {
    const report_section_t *g = r->sections + i;
    open_row(w);
    html_open(w, "td");
    html_open(w, "a");
    html_attr(w, "href", "#g%zu", i);
    html_text(w, g->name);
    html_close(w);
    html_close(w);
    td(w, "%zu", g->n_files);
    if (g->counts == NULL) {
      html_open(w, "td");
      html_attr(w, "colspan", "%d", r->references != NULL ? 4 : 3);
      html_text(w, "not processed yet");
      html_close(w);
    } else {
      const counter_reads_t *s = counter_reads(g->counts);
      td(w, "%s", html_num(s->reads_in).s);
      td(w, "%.1f%%", html_percent(s->reads_out, s->reads_in));
      td(w, "%.1f%%", html_percent(s->n_ends[2], s->reads_in));
      if (r->references != NULL)
        td(w, "%s", html_num(s->origin[COUNTER_LEFT][ORIGIN_INTERNAL] +
                             s->origin[COUNTER_RIGHT][ORIGIN_INTERNAL]).s);
    }
    html_close(w);
  }
  html_close(w);
}

/* the errors in a primer: the option, or decided automatically */
static void write_max_errors(html_writer_t *w, const report_t *r) {
  const report_setting_t *s = r->max_errors;
  double value = s != NULL ? s->value : r->options->max_errors;
  open_row(w);
  td(w, "errors in a primer");
  html_open(w, "td");
  html_attr(w, "class", "seq");
  if (s != NULL && s->pending)
    html_text(w, "automatic, at the first input file");
  else {
    if (value < 1)
      html_textf(w, "%g of its length", value);
    else
      html_textf(w, "%g", value);
    if (s != NULL && s->automatic) {
      html_text(w, " ");
      open_class(w, "span", "note");
      html_textf(w, "(automatic: error rate %.1f%% in the first %s reads "
                 "of %s)", 100 * s->measured, html_num(s->n_reads).s,
                 s->file);
      html_close(w);
    }
  }
  html_close(w);
  html_close(w);
}

/* the parameters of the run, at the right of the logo */
static void write_parameters(html_writer_t *w, const report_t *r) {
  const options_t *o = r->options;
  char n_refs[32], window[32], min_length[32], probe[32], probe_errors[32],
       terminal_window[32], threads[32];
  snprintf(n_refs, sizeof(n_refs), "%zu", r->references == NULL ? 0 :
           r->references->n_records);
  snprintf(window, sizeof(window), "%ld", o->window);
  snprintf(min_length, sizeof(min_length), "%zu", o->min_length);
  snprintf(probe, sizeof(probe), "%zu", o->probe_len);
  snprintf(probe_errors, sizeof(probe_errors), "%g", o->probe_error_rate);
  snprintf(terminal_window, sizeof(terminal_window), "%ld",
           o->terminal_window);
  snprintf(threads, sizeof(threads), "%zu", o->threads);
  html_open(w, "table");
  html_template(w, run_tmpl, "version", r->version,
                "command", r->command_line, "scheme", r->scheme->name,
                "policy", o->policy == TRIM_PRIMER ?
                  "primer (every primer-derived base removed)" :
                  "terminus (the conserved termini kept)",
                "references", r->references != NULL ? o->references : NULL,
                "n_references", n_refs,
                "no_references", r->references == NULL ? "" : NULL);
  write_max_errors(w, r);
  html_template(w, options_tmpl, "platform",
                o->illumina ? "illumina" : "nanopore", "window", window,
                "min_length", min_length,
                "probe", r->references != NULL ? probe : NULL,
                "probe_errors", probe_errors,
                "terminal_window", terminal_window,
                "search", o->exhaustive_search ? ", exhaustive search" : "",
                "threads", threads);
  html_close(w);
}

/* --- report --- */

static void write_report(html_writer_t *w, const report_t *r) {
  const html_logo_t logo = {report_logo, report_logo_dark, "flucut",
                            REPORT_LOGO_HEIGHT};
  html_page_start_logo(w, "flucut report", &logo, "report", r->refresh);
  write_parameters(w, r);
  html_page_head_end(w);
  html_template(w, style_tmpl);
  if (r->running) {
    open_class(w, "p", "running");
    html_textf(w, "Run in progress: the report is updated after each input "
               "file%s.", r->refresh > 0 ?
               " and this page is reloaded automatically" : "");
    html_close(w);
  }
  if (r->n_sections > 1) {
    html_tag_text(w, "h2", "Groups");
    write_groups_table(w, r);
  }
  html_tag_text(w, "h2", "All reads");
  write_section(w, r, r->total);
  if (r->n_sections > 1)
    for (size_t i = 0; i < r->n_sections; i++) {
      const report_section_t *g = r->sections + i;
      html_open(w, "h2");
      html_attr(w, "id", "g%zu", i);
      html_textf(w, "Group: %s", g->name);
      html_close(w);
      if (g->counts == NULL)
        note(w, "Not processed yet.");
      else
        write_section(w, r, g->counts);
    }
  html_page_end(w, "Generated by flucut "
                "(https://github.com/baia-ipc/flucut).");
}

/* --- dashboard (live mode): the same content on one screen --- */

/* the state and the parameters of the run, in the top bar of the
   dashboard (compact, one line each) */
static void write_dashboard_bar(html_writer_t *w, const report_t *r) {
  if (r->running) {
    open_class(w, "p", "running");
    html_text(w, "Run in progress: updated after each input file.");
    html_close(w);
  } else {
    html_open(w, "p");
    html_tag_text(w, "b", "End of the run.");
    html_close(w);
  }
  html_open(w, "p");
  html_textf(w, "%s \xc2\xb7 scheme %s \xc2\xb7 %s",
             r->options->illumina ? "Illumina" : "Nanopore", r->scheme->name,
             r->options->policy == TRIM_PRIMER ? "primer mode" :
             "terminus mode (conserved ends kept)");
  if (r->references != NULL)
    html_textf(w, " \xc2\xb7 %zu references", r->references->n_records);
  html_textf(w, " \xc2\xb7 flucut %s", r->version);
  html_close(w);
  const report_setting_t *e = r->max_errors;
  if (e != NULL && e->pending) {
    html_open(w, "p");
    html_text(w, "errors in a primer: automatic, at the first input file");
    html_close(w);
  } else if (e != NULL && e->automatic) {
    html_open(w, "p");
    html_textf(w, "errors in a primer: automatic, error rate %.1f%% in the "
               "first %s reads", 100 * e->measured,
               html_num(e->n_reads).s);
    html_close(w);
  }
}

/* a bar per barcode, its reads (of the run total) */
static void write_barcode_reads_bars(html_writer_t *w, const report_t *r) {
  size_t n = r->n_sections;
  const char **labels;
  size_t *values;
  const char **classes;
  XMALLOC(labels, sizeof(char *) * n);
  XMALLOC(values, sizeof(size_t) * n);
  XMALLOC(classes, sizeof(char *) * n);
  size_t total = 0;
  for (size_t i = 0; i < n; i++) {
    labels[i] = r->sections[i].name;
    values[i] = r->sections[i].counts != NULL ?
                counter_reads(r->sections[i].counts)->reads_in : 0;
    classes[i] = "s";
    total += values[i];
  }
  html_hbars(w, n, (const char *const *)labels, values,
             (const char *const *)classes, total, 0);
  free(labels);
  free(values);
  free(classes);
}

/* the panels of one sample (a barcode, or all the reads of a run of a
   single sample): the reads at each step, the internal priming (or, without
   references, the primers), the complete amplicons, and, across the whole
   width, the coverage of every segment by read origin. Each barcode is a
   different sample, so these are never pooled across barcodes */
static void write_dash_section(html_writer_t *w, const report_t *r,
    const counter_t *c) {
  const counter_reads_t *s = counter_reads(c);
  bool refs = r->references != NULL;
  open_class(w, "div", "tiles");
  write_tiles(w, r, s);
  html_close(w);
  write_stats(w, r, s);
  if (s->reads_in == 0) {
    note(w, "No reads yet.");
    return;
  }
  open_class(w, "div", "dgrid");
  html_panel_start(w, "Reads at each step", NULL);
  write_steps(w, s);
  html_panel_end(w);
  html_panel_start(w, "Primers found", NULL);
  write_primers(w, r, c);
  html_panel_end(w);
  if (refs) {
    html_panel_start(w, "Internal priming", NULL);
    write_sites(w, r, c);
    html_panel_end(w);
  }
  html_panel_start(w, "Complete amplicons", NULL);
  if (s->n_ends[2] > 0)
    write_lengths(w, r, c);
  else
    note(w, "None yet.");
  html_panel_end(w);
  if (refs) {
    html_panel_start(w, "Coverage by read origin", "full");
    write_coverage(w, r, c);
    html_panel_end(w);
  }
  html_close(w);   /* dgrid */
}

static void write_dashboard(html_writer_t *w, const report_t *r) {
  const html_logo_t logo = {report_logo, report_logo_dark, "flucut",
                            DASHBOARD_LOGO_HEIGHT};
  html_dashboard_start(w, "flucut dashboard", &logo, "dashboard",
                       r->refresh);
  write_dashboard_bar(w, r);
  html_dashboard_bar_end(w);
  html_template(w, style_tmpl);
  if (r->n_sections <= 1) {
    write_dash_section(w, r, r->total);
    html_page_end(w, NULL);
    return;
  }
  /* several barcode directories: each is a different sample, shown in its
     own tab (its own coverage and sites, not pooled with the others), with
     an overview tab of the run as a whole. The tabs are selected through
     the URL fragment, so the chosen one survives the automatic reload */
  open_class(w, "div", "tabs");
  html_open(w, "a");
  html_attr(w, "href", "#ov");
  html_text(w, "All reads");
  html_close(w);
  for (size_t i = 0; i < r->n_sections; i++) {
    html_open(w, "a");
    html_attr(w, "href", "#s%zu", i);
    html_text(w, r->sections[i].name);
    html_close(w);
  }
  html_close(w);
  /* the overview: the throughput of the whole run (which does not depend on
     the sample) and the per-barcode summary; no pooled coverage or sites */
  html_open(w, "section");
  html_attr(w, "id", "ov");
  html_attr(w, "class", "tabpanel");
  const counter_reads_t *s = counter_reads(r->total);
  open_class(w, "div", "tiles");
  write_tiles(w, r, s);
  html_close(w);
  write_stats(w, r, s);
  open_class(w, "div", "dgrid");
  html_panel_start(w, "Reads per barcode", NULL);
  if (s->reads_in > 0)
    write_barcode_reads_bars(w, r);
  else
    note(w, "No reads yet.");
  html_panel_end(w);
  html_panel_start(w, "Reads at each step (all barcodes)", "wide");
  if (s->reads_in > 0)
    write_steps(w, s);
  else
    note(w, "No reads yet.");
  html_panel_end(w);
  html_panel_start(w, "Barcodes", "full");
  write_groups_table(w, r);
  html_panel_end(w);
  html_close(w);   /* dgrid */
  html_close(w);   /* section#ov */
  for (size_t i = 0; i < r->n_sections; i++) {
    const report_section_t *g = r->sections + i;
    html_open(w, "section");
    html_attr(w, "id", "s%zu", i);
    html_attr(w, "class", "tabpanel bc");
    if (g->counts == NULL)
      note(w, "Not processed yet.");
    else
      write_dash_section(w, r, g->counts);
    html_close(w);
  }
  html_page_end(w, NULL);
}

void html_report_write(const char *path, const report_t *report,
    Logger *logger) {
  assert(path != NULL);
  assert(report != NULL);
  assert(report->total != NULL);
  FILE *f = xfopen(path, REPORT_TMP_SFX, "w", "HTML report", logger);
  html_writer_t w;
  html_writer_init(&w, f);
  if (report->dashboard)
    write_dashboard(&w, report);
  else
    write_report(&w, report);
  xfclose(&f, logger);
  xfrename(path, REPORT_TMP_SFX, path, NULL, logger);
}
