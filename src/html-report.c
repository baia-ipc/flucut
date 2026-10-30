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

static void write_section(html_writer_t *w, const report_t *r,
    const counter_t *c) {
  const counter_reads_t *s = counter_reads(c);
  open_class(w, "div", "tiles");
  write_tiles(w, r, s);
  html_close(w);
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

void html_report_write(const char *path, const report_t *report,
    Logger *logger) {
  assert(path != NULL);
  assert(report != NULL);
  assert(report->total != NULL);
  FILE *f = xfopen(path, REPORT_TMP_SFX, "w", "HTML report", logger);
  html_writer_t w;
  html_writer_init(&w, f);
  write_report(&w, report);
  xfclose(&f, logger);
  xfrename(path, REPORT_TMP_SFX, path, NULL, logger);
}
