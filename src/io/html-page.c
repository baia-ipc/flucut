/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "html-page.h"
#include "html-template.h"
#include "html-writer.h"
#include "numutils.h"
#include "strbuf.h"
#include "test.h"

/* plots: width (the height depends on the plot) */
#define HTML_PLOT_WIDTH 760
/* plots with axes: margins of the plot area (the left one for the labels
   of the y axis) */
#define HTML_PLOT_LEFT 60
#define HTML_PLOT_RIGHT 20
#define HTML_PLOT_TOP 10
/* plots with axes: distance from the x axis of the tick labels */
#define HTML_TICK_LABEL_DY 16
/* bar chart: width of the labels, space right of the bars for the values,
   height of the bars and space between them */
#define HTML_LABEL_WIDTH 190
#define HTML_VALUE_WIDTH 150
#define HTML_BAR_HEIGHT 16
#define HTML_BAR_GAP 4
/* histogram: height of the plot area, space below it (x axis labels and
   title), distance of the title from the x axis */
#define HTML_HIST_HEIGHT 200
#define HTML_HIST_BOTTOM 30
#define HTML_HIST_TITLE_DY 28
/* rank plot: the same; maximal number of points of the line */
#define HTML_KNEE_HEIGHT 240
#define HTML_KNEE_BOTTOM 34
#define HTML_KNEE_TITLE_DY 30
#define HTML_MAX_KNEE_POINTS 400

/* the style and the templates of the pages, embedded by the Makefile */
const char *const html_page_css =
#include "html-page-css.inc"
  ;
static const char *const start_tmpl =
#include "html-page-start-html.inc"
  ;
static const char *const end_tmpl =
#include "html-page-end-html.inc"
  ;
static const char *const brand_tmpl =
#include "html-page-brand-html.inc"
  ;
static const char *const dashboard_brand_tmpl =
#include "html-page-dashboard-brand-html.inc"
  ;

/* the time of writing, in <date> (at least 32 bytes); empty if not
   available; localtime_r: the pages can be written by several threads */
static void now_date(char *date, size_t size) {
  time_t now = time(NULL);
  struct tm now_tm;
  if (localtime_r(&now, &now_tm) == NULL ||
      strftime(date, size, "%Y-%m-%d %H:%M:%S", &now_tm) == 0)
    date[0] = '\0';
}

/* the head of the page and the start of the body (<updated>: the time of
   writing, for the reload) */
static void write_start(html_writer_t *w, const char *title, long refresh,
    const char *updated, const char *body_class) {
  assert(w->depth == 0);
  char max[32], min[32];
  snprintf(max, sizeof(max), "%ld", refresh);
  snprintf(min, sizeof(min), "%d", HTML_PAGE_MIN_RELOAD);
  html_template(w, start_tmpl, "title", title,
                "refresh", refresh > 0 ? max : NULL, "min_refresh", min,
                "updated", updated, "css", html_page_css,
                "body_class", body_class);
}

/* the logo (NULL: none) and the heading over the time of writing */
static void write_brand(html_writer_t *w, const char *tmpl,
    const html_logo_t *logo, const char *heading, const char *updated) {
  assert(logo == NULL || logo->src != NULL);
  char height[32];
  snprintf(height, sizeof(height), "%d", logo != NULL ? logo->height : 0);
  html_template(w, tmpl,
                "logo", logo != NULL ? logo->src : NULL,
                "logo_dark", logo != NULL ? logo->src_dark : NULL,
                "logo_alt", logo == NULL ? NULL :
                  logo->alt != NULL ? logo->alt : "",
                "logo_height", logo != NULL ? height : NULL,
                "heading", heading, "updated", updated);
}

void html_page_start(html_writer_t *w, const char *title, long refresh) {
  assert(w != NULL);
  assert(title != NULL);
  char date[64];
  now_date(date, sizeof(date));
  write_start(w, title, refresh, date, "page");
  html_tag_text(w, "h1", title);
  html_open(w, "p");
  html_attr(w, "class", "sub updated");
  html_textf(w, "Updated: %s", date);
  html_close(w);
}

void html_page_start_logo(html_writer_t *w, const char *title,
    const html_logo_t *logo, const char *heading, long refresh) {
  assert(w != NULL);
  assert(title != NULL);
  assert(logo != NULL);
  assert(heading != NULL);
  char date[64];
  now_date(date, sizeof(date));
  write_start(w, title, refresh, date, "page");
  html_open(w, "div");
  html_attr(w, "class", "head");
  write_brand(w, brand_tmpl, logo, heading, date);
  html_open(w, "div");
  html_attr(w, "class", "side");
}

void html_page_head_end(html_writer_t *w) {
  assert(w != NULL);
  html_close(w);  /* div.side */
  html_close(w);  /* div.head */
}

void html_dashboard_start(html_writer_t *w, const char *title,
    const html_logo_t *logo, const char *heading, long refresh) {
  assert(w != NULL);
  assert(title != NULL);
  assert(heading != NULL);
  char date[64];
  now_date(date, sizeof(date));
  write_start(w, title, refresh, date, "dash");
  html_open(w, "header");
  html_attr(w, "class", "dbar");
  write_brand(w, dashboard_brand_tmpl, logo, heading, date);
  html_open(w, "div");
  html_attr(w, "class", "dinfo");
}

void html_dashboard_bar_end(html_writer_t *w) {
  assert(w != NULL);
  html_close(w);  /* div.dinfo */
  html_close(w);  /* header.dbar */
}

void html_panel_start(html_writer_t *w, const char *title,
    const char *cls) {
  assert(w != NULL);
  assert(title != NULL);
  html_open(w, "section");
  html_attr(w, "class", "panel%s%s", cls != NULL ? " " : "",
            cls != NULL ? cls : "");
  html_tag_text(w, "h3", title);
}

void html_panel_end(html_writer_t *w) {
  assert(w != NULL);
  html_close(w);
}

void html_page_end(html_writer_t *w, const char *footer) {
  assert(w != NULL);
  assert(w->depth == 0);
  html_template(w, end_tmpl, "footer", footer);
  html_writer_end(w);
}

void html_tile(html_writer_t *w, const char *value, const char *label) {
  assert(w != NULL);
  assert(value != NULL);
  assert(label != NULL);
  html_open(w, "div");
  html_attr(w, "class", "tile");
  html_open(w, "div");
  html_attr(w, "class", "v");
  html_text(w, value);
  html_close(w);
  html_open(w, "div");
  html_attr(w, "class", "l");
  html_text(w, label);
  html_close(w);
  html_close(w);
}

html_num_t html_num(size_t n) {
  html_num_t out;
  char digits[32];
  int len = snprintf(digits, sizeof(digits), "%zu", n), j = 0;
  /* (at most 20 digits and 6 separators) */
  for (int i = 0; i < len; i++) {
    if (i > 0 && (len - i) % 3 == 0)
      out.s[j++] = ',';
    out.s[j++] = digits[i];
  }
  out.s[j] = '\0';
  return out;
}

double html_percent(size_t a, size_t b) {
  return b == 0 ? 0.0 : 100.0 * (double)a / (double)b;
}

/* --- plots --- */

/* start of a plot of the given height (the width is HTML_PLOT_WIDTH) */
static void svg_start(html_writer_t *w, int height) {
  html_open(w, "svg");
  html_attr(w, "viewBox", "0 0 %d %d", HTML_PLOT_WIDTH, height);
  html_attr(w, "width", "%d", HTML_PLOT_WIDTH);
  html_attr(w, "role", "img");
}

/* a text label at (<x>, <y>) of the class <cls> (NULL: none) aligned as
   <anchor> (NULL: start), with the text (formatted as by printf) */
static void svg_text(html_writer_t *w, const char *cls, double x, double y,
    const char *anchor, const char *fmt, ...)
  __attribute__((format(printf, 6, 7)));

static void svg_text(html_writer_t *w, const char *cls, double x, double y,
    const char *anchor, const char *fmt, ...) {
  html_open(w, "text");
  if (cls != NULL)
    html_attr(w, "class", "%s", cls);
  html_attr(w, "x", "%.1f", x);
  html_attr(w, "y", "%.1f", y);
  if (anchor != NULL)
    html_attr(w, "text-anchor", "%s", anchor);
  va_list ap;
  va_start(ap, fmt);
  html_vtextf(w, fmt, ap);
  va_end(ap);
  html_close(w);
}

/* a line from (<x1>, <y1>) to (<x2>, <y2>) of the class <cls> (more
   attributes can follow) */
static void svg_line(html_writer_t *w, const char *cls, double x1,
    double y1, double x2, double y2) {
  html_elem(w, "line");
  html_attr(w, "class", "%s", cls);
  html_attr(w, "x1", "%.1f", x1);
  html_attr(w, "y1", "%.1f", y1);
  html_attr(w, "x2", "%.1f", x2);
  html_attr(w, "y2", "%.1f", y2);
}

void html_hbars(html_writer_t *w, size_t n, const char *const *labels,
    const size_t *values, const char *const *classes, size_t total,
    double mean) {
  assert(w != NULL);
  assert(n == 0 || (labels != NULL && values != NULL && classes != NULL));
  size_t max = 1;
  for (size_t i = 0; i < n; i++)
    if (values[i] > max)
      max = values[i];
  int area = HTML_PLOT_WIDTH - HTML_LABEL_WIDTH - HTML_VALUE_WIDTH;
  int height = (int)n * (HTML_BAR_HEIGHT + HTML_BAR_GAP) + 10;
  svg_start(w, height);
  if (mean > 0) {
    double x = HTML_LABEL_WIDTH + (double)area * mean / (double)max;
    html_open(w, "line");
    html_attr(w, "class", "a");
    html_attr(w, "x1", "%.1f", x);
    html_attr(w, "y1", "0");
    html_attr(w, "x2", "%.1f", x);
    html_attr(w, "y2", "%d", height);
    html_attr(w, "stroke-dasharray", "4 3");
    html_inline(w);
    html_tag_text(w, "title", "mean");
    html_close(w);
  }
  for (size_t i = 0; i < n; i++) {
    int y = (int)i * (HTML_BAR_HEIGHT + HTML_BAR_GAP) + 4,
        text_y = y + HTML_BAR_HEIGHT - 4;
    double bar = (double)area * (double)values[i] / (double)max;
    svg_text(w, NULL, HTML_LABEL_WIDTH - 8, text_y, "end", "%s", labels[i]);
    /* (a bar of a value > 0 is visible) */
    html_open(w, "rect");
    html_attr(w, "class", "%s", classes[i]);
    html_attr(w, "x", "%d", HTML_LABEL_WIDTH);
    html_attr(w, "y", "%d", y);
    html_attr(w, "width", "%.1f", bar > 0 && bar < 1 ? 1.0 : bar);
    html_attr(w, "height", "%d", HTML_BAR_HEIGHT);
    html_attr(w, "rx", "2");
    html_inline(w);
    html_tag_textf(w, "title", "%s: %s", labels[i], html_num(values[i]).s);
    html_close(w);
    svg_text(w, "muted halo", HTML_LABEL_WIDTH + bar + 5, text_y, NULL,
             "%s (%.1f%%)", html_num(values[i]).s,
             html_percent(values[i], total));
  }
  html_close(w);
}

void html_histogram(html_writer_t *w, const size_t *hist, size_t n_bins_all,
    size_t bin_width, const char *what, const char *unit,
    const char *empty) {
  assert(w != NULL);
  assert(hist != NULL || n_bins_all == 0);
  assert(bin_width > 0);
  assert(what != NULL);
  assert(unit != NULL);
  assert(empty != NULL);
  size_t total = 0, max = 0, last = 0;
  for (size_t b = 0; b < n_bins_all; b++) {
    total += hist[b];
    if (hist[b] > max)
      max = hist[b];
    if (hist[b] > 0)
      last = b;
  }
  if (total == 0) {
    html_open(w, "p");
    html_attr(w, "class", "note");
    html_text(w, empty);
    html_close(w);
    return;
  }
  /* show the bins up to the 99.5th percentile */
  size_t cum = 0, shown = 0;
  for (size_t b = 0; b <= last; b++) {
    cum += hist[b];
    shown = b;
    if ((double)cum >= 0.995 * (double)total)
      break;
  }
  size_t n_bins = shown + 1;
  int left = HTML_PLOT_LEFT, bottom = HTML_HIST_BOTTOM, top = HTML_PLOT_TOP,
      width = HTML_PLOT_WIDTH - left - HTML_PLOT_RIGHT,
      height = HTML_HIST_HEIGHT, x_axis = top + height;
  double ymax = nice_number((double)max), bw = (double)width / (double)n_bins;
  svg_start(w, height + top + bottom);
  /* horizontal grid, with the counts */
  for (int g = 0; g <= 4; g++) {
    double y = x_axis - height * g / 4.0;
    svg_line(w, "grid", left, y, left + width, y);
    svg_text(w, "muted", left - 6, y + 4, "end", "%s",
             html_num((size_t)(ymax * g / 4)).s);
  }
  for (size_t b = 0; b < n_bins; b++) {
    double h = height * (double)hist[b] / ymax;
    html_open(w, "rect");
    html_attr(w, "class", "s");
    html_attr(w, "x", "%.1f", left + b * bw);
    html_attr(w, "y", "%.1f", x_axis - h);
    html_attr(w, "width", "%.1f", bw > 2 ? bw - 1 : bw);
    html_attr(w, "height", "%.1f", h);
    html_inline(w);
    /* (the last bin: all larger values) */
    if (b == n_bins_all - 1)
      html_tag_textf(w, "title", "%zu %s or more: %s", b * bin_width, unit,
                     html_num(hist[b]).s);
    else
      html_tag_textf(w, "title", "%zu-%zu %s: %s", b * bin_width,
                     (b + 1) * bin_width - 1, unit, html_num(hist[b]).s);
    html_close(w);
  }
  /* x axis: ticks and values */
  size_t max_len = n_bins * bin_width,
         step = (size_t)nice_number((double)max_len / 8);
  for (size_t x = 0; x <= max_len; x += step) {
    double px = left + (double)x / bin_width * bw;
    svg_line(w, "axis", px, x_axis, px, x_axis + 4);
    svg_text(w, "muted", px, x_axis + HTML_TICK_LABEL_DY, "middle", "%zu",
             x);
  }
  svg_text(w, "muted", left + width / 2, x_axis + HTML_HIST_TITLE_DY,
           "middle", "%s (%s)%s", what, unit,
           shown < last ? ", longest 0.5% not shown" : "");
  html_close(w);
}

void html_knee(html_writer_t *w, const uint64_t *counts, size_t n,
    size_t n_selected, const char *x_label, const char *y_label,
    const char *selected_label) {
  assert(w != NULL);
  assert(counts != NULL || n == 0);
  assert(x_label != NULL);
  assert(y_label != NULL);
  assert(selected_label != NULL);
  if (n == 0 || counts[0] == 0)
    return;
  int left = HTML_PLOT_LEFT, bottom = HTML_KNEE_BOTTOM, top = HTML_PLOT_TOP,
      width = HTML_PLOT_WIDTH - left - HTML_PLOT_RIGHT,
      height = HTML_KNEE_HEIGHT, x_axis = top + height;
  double lx = log10((double)n) > 0 ? log10((double)n) : 1,
         ly = log10((double)counts[0]) > 0 ? log10((double)counts[0]) : 1;
  svg_start(w, height + top + bottom);
  /* grid at the powers of 10 of the counts, with their values */
  for (int e = 0; e <= (int)ceil(ly); e++) {
    double y = x_axis - height * e / ly;
    if (y < top)
      break;
    svg_line(w, "grid", left, y, left + width, y);
    svg_text(w, "muted", left - 6, y + 4, "end", "%s",
             html_num((size_t)pow(10, e)).s);
  }
  /* the ranks at the powers of 10 */
  for (int e = 0; e <= (int)ceil(lx); e++) {
    double x = left + width * e / lx;
    if (x > left + width + 1)
      break;
    svg_text(w, "muted", x, x_axis + HTML_TICK_LABEL_DY, "middle", "%s",
             html_num((size_t)pow(10, e)).s);
  }
  /* the line, with at most about HTML_MAX_KNEE_POINTS points */
  strbuf_t points;
  strbuf_init(&points);
  double prev = -1;
  for (size_t i = 0; i < n; i++) {
    double pos = log10((double)(i + 1)) / lx;
    if (i + 1 < n && pos - prev < 1.0 / HTML_MAX_KNEE_POINTS)
      continue;
    prev = pos;
    double y = counts[i] > 0 ? log10((double)counts[i]) / ly : 0;
    char point[64];
    snprintf(point, sizeof(point), "%s%.1f,%.1f", points.len > 0 ? " " : "",
             left + width * pos, x_axis - height * y);
    strbuf_append_str(&points, point);
  }
  html_elem(w, "polyline");
  html_attr(w, "class", "line");
  html_attr(w, "points", "%.*s", (int)points.len, points.s);
  strbuf_free(&points);
  if (n_selected > 0) {
    /* the line at most at the last rank, the label with the number */
    size_t rank = n_selected < n ? n_selected : n;
    double x = left + width * log10((double)rank) / lx;
    svg_line(w, "a", x, top, x, x_axis);
    html_attr(w, "stroke-dasharray", "4 3");
    svg_text(w, NULL, x + 4, top + 12, NULL, "%s %s",
             html_num(n_selected).s, selected_label);
  }
  svg_text(w, "muted", left + width / 2, x_axis + HTML_KNEE_TITLE_DY,
           "middle", "%s (log scale)", x_label);
  html_open(w, "text");
  html_attr(w, "class", "muted");
  html_attr(w, "transform", "translate(14,%d) rotate(-90)", top + height / 2);
  html_attr(w, "text-anchor", "middle");
  html_textf(w, "%s (log scale)", y_label);
  html_close(w);
  html_close(w);
}

/* --- unit tests --- */

void html_page_test(void) {
  test_eq_str(html_num(0).s, "0");
  test_eq_str(html_num(999).s, "999");
  test_eq_str(html_num(1234567).s, "1,234,567");
  test_eq_str(html_num(SIZE_MAX).s, SIZE_MAX == UINT64_MAX ?
              "18,446,744,073,709,551,615" : html_num(SIZE_MAX).s);
  test_assert(html_percent(1, 4) == 25.0 && html_percent(1, 0) == 0.0);
  char *out = NULL;
  size_t size = 0;
  FILE *f = open_memstream(&out, &size);
  test_assert(f != NULL);
  html_writer_t w;
  html_writer_init(&w, f);
  html_page_start(&w, "Test <report>", 30);
  html_open(&w, "div");
  html_attr(&w, "class", "tiles");
  html_tile(&w, "1,234", "reads & more");
  html_close(&w);
  const char *labels[] = {"a", "other"};
  const size_t values[] = {30, 10};
  const char *classes[] = {"s", "p"};
  html_hbars(&w, 2, labels, values, classes, 40, 20);
  const size_t hist[] = {0, 2, 5, 1};
  html_histogram(&w, hist, 4, 10, "length", "bp", "No values.");
  const size_t none[] = {0, 0};
  html_histogram(&w, none, 2, 10, "length", "bp", "No values.");
  const uint64_t counts[] = {500, 300, 100, 5, 1};
  html_knee(&w, counts, 5, 3, "items, by decreasing count", "count",
            "selected");
  /* more selected than counts: the line at the last rank */
  html_knee(&w, counts, 5, 9, "items", "count", "selected");
  html_page_end(&w, "Generated by test");
  fclose(f);
  const char *exp[] = {"<!doctype html>\n<html lang=\"en\">\n<head>\n",
    "  <title>Test &lt;report&gt;</title>\n",
    "  <noscript><meta http-equiv=\"refresh\" content=\"30\"></noscript>\n",
    "      var max = 30, min = 1, updated = \"", "sessionStorage",
    "  <style>\n    :root{", "\n  </style>\n</head>\n<body class=\"page\">\n"
    "<main>\n<h1>Test &lt;report&gt;</h1>\n<p class=\"sub updated\">Updated: ",
    "<div class=\"tiles\">\n  <div class=\"tile\">\n"
    "    <div class=\"v\">1,234</div>\n"
    "    <div class=\"l\">reads &amp; more</div>\n  </div>\n</div>\n",
    "<svg viewBox=\"0 0 760 50\" width=\"760\" role=\"img\">\n",
    "  <line class=\"a\" x1=\"470.0\" y1=\"0\" x2=\"470.0\" y2=\"50\" "
    "stroke-dasharray=\"4 3\"><title>mean</title></line>\n",
    "  <text x=\"182.0\" y=\"16.0\" text-anchor=\"end\">a</text>\n"
    "  <rect class=\"s\" x=\"190\" y=\"4\" width=\"420.0\" height=\"16\" "
    "rx=\"2\"><title>a: 30</title></rect>\n"
    "  <text class=\"muted halo\" x=\"615.0\" y=\"16.0\">30 (75.0%)</text>\n",
    "<title>20-29 bp: 5</title>", "<title>30 bp or more: 1</title>",
    "  <line class=\"grid\" x1=\"60.0\" y1=\"210.0\" x2=\"740.0\" "
    "y2=\"210.0\"/>\n",
    ">length (bp)</text>\n</svg>\n",
    "<p class=\"note\">No values.</p>\n",
    "  <polyline class=\"line\" points=\"60.0,10.0 ", "3 selected</text>",
    "  <line class=\"a\" x1=\"740.0\" y1=\"10.0\"", "9 selected</text>",
    ">items, by decreasing count (log scale)</text>",
    ">count (log scale)</text>\n</svg>\n",
    "</svg>\n<p class=\"note\">Generated by test</p>\n</main>\n</body>\n"
    "</html>\n"};
  for (size_t i = 0; i < sizeof(exp) / sizeof(exp[0]); i++)
    if (strstr(out, exp[i]) == NULL) {
      fprintf(stderr, "Not found in the page: \"%s\"\n", exp[i]);
      test_assert(false);
    }
  free(out);
  /* with a logo */
  f = open_memstream(&out, &size);
  test_assert(f != NULL);
  html_writer_init(&w, f);
  html_logo_t logo = {"logo.png", "logo-dark.png", "Test", 80};
  html_page_start_logo(&w, "Test report", &logo, "report", 0);
  html_tag_text(&w, "p", "side");
  html_page_head_end(&w);
  html_page_end(&w, NULL);
  fclose(f);
  test_assert(strstr(out, "<title>Test report</title>") != NULL);
  test_assert(strstr(out, "http-equiv") == NULL);
  test_assert(strstr(out, "<div class=\"head\">\n  <div class=\"brand\">\n"
              "    <h1 class=\"logo\">\n      <picture>\n        <source "
              "media=\"(prefers-color-scheme: dark)\" srcset=\""
              "logo-dark.png\">\n        <img src=\"logo.png\" alt=\"Test\" "
              "height=\"80\">\n      </picture>\n      <span>report</span>\n"
              "    </h1>\n    <p class=\"sub updated\">Updated: ") != NULL);
  test_assert(strstr(out, "  <div class=\"side\">\n    <p>side</p>\n"
              "  </div>\n</div>\n</main>\n") != NULL);
  free(out);
  /* dashboard, without a logo */
  f = open_memstream(&out, &size);
  test_assert(f != NULL);
  html_writer_init(&w, f);
  html_dashboard_start(&w, "Test dashboard", NULL, "dashboard", 30);
  html_tag_text(&w, "p", "state");
  html_dashboard_bar_end(&w);
  html_open(&w, "div");
  html_attr(&w, "class", "dgrid");
  html_panel_start(&w, "Values", "tall");
  html_panel_end(&w);
  html_panel_start(&w, "Steps", NULL);
  html_tag_text(&w, "p", "x");
  html_panel_end(&w);
  html_close(&w);
  html_page_end(&w, NULL);
  fclose(f);
  test_assert(strstr(out, "<body class=\"dash\">") != NULL);
  test_assert(strstr(out, "<picture>") == NULL);
  test_assert(strstr(out, "<header class=\"dbar\">\n  <div class=\"brand\">\n"
              "    <span>dashboard</span>\n") != NULL);
  test_assert(strstr(out, "  <div class=\"dinfo\">\n    <p>state</p>\n"
              "  </div>\n</header>\n") != NULL);
  test_assert(strstr(out, "<div class=\"dgrid\">\n"
              "  <section class=\"panel tall\">\n    <h3>Values</h3>\n"
              "  </section>\n  <section class=\"panel\">\n    <h3>Steps</h3>\n"
              "    <p>x</p>\n  </section>\n</div>\n") != NULL);
  free(out);
}
