/* (c) 2026 Giorgio Gonnella */

/*
   Self-contained HTML reports: a page with a style for light and dark
   themes (no external resources), tiles with the main values, and plots as
   inline SVG (horizontal bar chart, histogram, rank plot).

   A report is written by a html_writer_t (html-writer.h) as:
   html_page_start(), the content (headings, paragraphs, tables, tiles,
   plots), html_page_end(). The style defines:
   - text: h1, h2, h3; p.sub (subtitle), p.note (small, muted), p.running
     (highlighted box); tables (numbers right-aligned, td.seq: monospace);
   - div.tiles containing the tiles (html_tile());
   - SVG classes: s (main color, e.g. samples), p (secondary, muted color),
     a (accent line), line, grid, axis, muted, halo; the colors are the
     CSS variables --primary (s, line), --secondary (p) and --accent (a).

   The style is in html-page.css, the static parts of the pages (head,
   heading with the logo, end) are templates (html-template.h) in
   html-page-*.html: the Makefile embeds them in html-page.c.
*/

#ifndef HTML_PAGE_H
#define HTML_PAGE_H

#include <stddef.h>
#include <stdint.h>
#include "html-writer.h"

/* the CSS of the page (written by html_page_start()) */
extern const char *const html_page_css;

/* minimum interval of the reload of a page (seconds) */
#define HTML_PAGE_MIN_RELOAD 1

/*
   Start the page (the writer must be at the top level): head (title,
   style; if <refresh> > 0, the page is reloaded, e.g. for a report
   rewritten during a run: at an interval halved, down to
   HTML_PAGE_MIN_RELOAD seconds, if the page has changed since the last
   reload, otherwise doubled, up to <refresh> seconds; this requires
   JavaScript, otherwise every <refresh> seconds), the title as h1 heading
   and the time of writing ("Updated: ...")
*/
void html_page_start(html_writer_t *w, const char *title, long refresh);

/* A logo: image URL (e.g. a data URI, "data:image/webp;base64,..."), an
   optional variant for the dark color scheme (NULL: none), the alternative
   text and the height of the image (pixels) */
typedef struct {
  const char *src, *src_dark, *alt;
  int height;
} html_logo_t;

/*
   As html_page_start(), with a head made of two columns: on the left, the
   h1 heading, made of the <logo> with <heading> under it (e.g. the title
   without the name shown by the logo; bold spaced capitals), centered
   over the time of writing; on the right (vertically centered), the
   content written next (e.g. a table of the parameters), up to
   html_page_head_end(); on narrow screens, the right column is moved
   below the left one
*/
void html_page_start_logo(html_writer_t *w, const char *title,
    const html_logo_t *logo, const char *heading, long refresh);

/* End the head started by html_page_start_logo() */
void html_page_head_end(html_writer_t *w);

/*
   Start a dashboard: a page made to be seen on a screen, e.g. during a
   run: it fits the window (on a large screen, the page does not scroll:
   the panels do, if needed), with a top bar: the <logo> (NULL: none) with
   <heading> under it and the time of writing, then the content written
   next (e.g. the state and parameters of the run), up to
   html_dashboard_bar_end(); then e.g. tiles (div.tiles) and a grid of
   panels (div.dgrid, see html_panel_start()); <refresh>: see
   html_page_start(); end with html_page_end()
*/
void html_dashboard_start(html_writer_t *w, const char *title,
    const html_logo_t *logo, const char *heading, long refresh);

/* End the top bar of a dashboard */
void html_dashboard_bar_end(html_writer_t *w);

/*
   Start a panel of a dashboard grid (in a div.dgrid: three columns, rows
   of equal height), with a <title>; <cls> (NULL: none): "tall" (two rows)
   or "wide" (two columns); the content scrolls if needed
*/
void html_panel_start(html_writer_t *w, const char *title, const char *cls);

/* End a panel */
void html_panel_end(html_writer_t *w);

/* End the page, with a <footer> note (NULL: none); all elements opened
   by the content must be closed (it ends the writer) */
void html_page_end(html_writer_t *w, const char *footer);

/* A tile with a value and its label (text), in a div.tiles */
void html_tile(html_writer_t *w, const char *value, const char *label);

/* a number formatted by html_num() */
typedef struct {
  char s[32];
} html_num_t;

/*
   A number with thousands separators (e.g. 1,234,567): the string is
   html_num(n).s, valid up to the end of the statement (a value, not a
   static buffer: any number of calls can be used in one statement, e.g.
   html_textf(w, "%s of %s", html_num(a).s, html_num(b).s))
*/
html_num_t html_num(size_t n);

/* 100 * a / b (0 if b is 0) */
double html_percent(size_t a, size_t b);

/*
   Horizontal bar chart: for each of the <n> values, its label, a bar of
   the class classes[i] (e.g. "s" or "p") and the value with the percentage
   of <total>; if <mean> > 0, a dashed line at this value
*/
void html_hbars(html_writer_t *w, size_t n, const char *const *labels,
    const size_t *values, const char *const *classes, size_t total,
    double mean);

/*
   Histogram of the counts of <n_bins> bins of width <bin_width> (> 0; bin
   b: values from b * bin_width to (b + 1) * bin_width - 1; the last bin:
   all larger values), up to the bin of the 99.5th percentile; x axis:
   "<what> (<unit>)"; if all counts are 0, the note <empty> instead
*/
void html_histogram(html_writer_t *w, const size_t *hist, size_t n_bins,
    size_t bin_width, const char *what, const char *unit,
    const char *empty);

/*
   Rank plot (log-log, "knee plot"): the <n> counts of the items, in
   decreasing order (e.g. the reads of each barcode); axes titles:
   "<x_label> (log scale)" and "<y_label> (log scale)"; if <n_selected> >
   0, a dashed line at that rank (at the last rank, <n>, if larger),
   labeled "<n_selected> <selected_label>" (e.g. the barcodes identified as
   cells: "cells"); the labels are text (HTML-escaped)
*/
void html_knee(html_writer_t *w, const uint64_t *counts, size_t n,
    size_t n_selected, const char *x_label, const char *y_label,
    const char *selected_label);

/* Unit tests */
void html_page_test(void);

#endif
