/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "alloc.h"
#include "html-page.h"
#include "html-report-markup.h"
#include "html-report-plots.h"
#include "svg-writer.h"
#include "segments.h"
#include "strbuf.h"

/* --- coverage by read origin --- */

/* the bands of the plot: the classes of the reads (COVER_*) gathered in
   the categories of the output (the strands together), in the order in
   which they are stacked, from the bottom: the reads which cover the
   whole segment evenly first (fragments, full templates), so that
   their bands are flat and the step of the internally primed reads above
   them is seen at its site; the reads with an internal deletion, which
   cover the two ends only, on top */
enum { BAND_FRAGMENT, BAND_FULL, BAND_INTERNAL, BAND_DELETION, N_BANDS };
static const int band_of_cover[COVER_N_CATEGORIES] = {
  BAND_FULL, BAND_FULL, BAND_INTERNAL, BAND_INTERNAL, BAND_FRAGMENT,
  BAND_FRAGMENT, BAND_DELETION};
static const char *const band_class[N_BANDS] = {
  "c-frag", "c-full", "c-int", "c-del"};
static const char *const band_label[N_BANDS] = {
  "fragment", "full template", "internally primed", "internal deletion"};

void plots_coverage_legend(html_writer_t *w, size_t min_depth) {
  char low[40];
  snprintf(low, sizeof(low), "below %zu\xc3\x97", min_depth);
  markup_open(w, "div", "lg");
  for (int b = BAND_FULL; b < N_BANDS; b++)
    markup_swatch(w, band_class[b], band_label[b]);
  markup_swatch(w, band_class[BAND_FRAGMENT], band_label[BAND_FRAGMENT]);
  markup_swatch(w, "c-low", low);
  html_close(w);
}

/* geometry of the coverage plot (pixels of the viewBox) */
enum { COV_W = 560, COV_H = 122, COV_LEFT = 52, COV_RIGHT = 8, COV_TOP = 12,
       COV_PLOT_H = 92, COV_MAX_POINTS = 260 };

/* append to <pts> the points of the profile <y> at the <n> abscissas <x>,
   forward or backward; a point between two of the same height is left
   out (the profiles are made of steps: most points are) */
static void append_profile(strbuf_t *pts, const double *x, const double *y,
    size_t n, bool backward) {
  char buf[40];
  for (size_t k = 0; k < n; k++) {
    size_t i = backward ? n - 1 - k : k;
    if (i > 0 && i < n - 1 && y[i] == y[i - 1] && y[i] == y[i + 1])
      continue;
    int len = snprintf(buf, sizeof(buf), "%s%.1f,%.1f",
                       pts->len > 0 ? " " : "", x[i], y[i]);
    strbuf_append(pts, buf, (size_t)len);
  }
}

/* an internal priming site is marked from this number of read ends */
#define MIN_SITE_READS 10

/* the position (0-based) of the reference <ref>, of length <len>, at
   which most internally primed read ends start; false if fewer than
   MIN_SITE_READS start at any */
static bool main_site(const counter_t *c, size_t ref, size_t len,
    size_t *pos) {
  const size_t *left = counter_sites(c, ref, COUNTER_LEFT);
  const size_t *right = counter_sites(c, ref, COUNTER_RIGHT);
  size_t best = 0;
  for (size_t p = 0; p < len; p++) {
    size_t n = (left != NULL ? left[p] : 0) + (right != NULL ? right[p] : 0);
    if (n > best) {
      best = n;
      *pos = p;
    }
  }
  return best >= MIN_SITE_READS;
}

/* the stacked bands with the depths <depth> (by band and position),
   scaled to the greatest total depth <max> */
static void coverage_bands(html_writer_t *w, const size_t *depth,
    size_t len, size_t max) {
  size_t n = len < COV_MAX_POINTS ? len : COV_MAX_POINTS;
  double inner = COV_W - COV_LEFT - COV_RIGHT;
  double *x, *low, *high;
  XMALLOC(x, sizeof(double) * n);
  XMALLOC(low, sizeof(double) * n);
  XMALLOC(high, sizeof(double) * n);
  for (size_t i = 0; i < n; i++) {
    x[i] = COV_LEFT + (double)i * inner / (double)(n - 1);
    low[i] = COV_TOP + COV_PLOT_H;
  }
  strbuf_t pts;
  strbuf_init(&pts);
  for (size_t band = 0; band < N_BANDS; band++) {
    const size_t *d = depth + band * len;
    bool any = false;
    for (size_t i = 0; i < n; i++) {
      size_t p = (size_t)((double)i * (double)(len - 1) / (double)(n - 1));
      high[i] = low[i] - COV_PLOT_H * (double)d[p] / (double)max;
      any = any || d[p] > 0;
    }
    if (!any)
      continue;
    strbuf_clear(&pts);
    append_profile(&pts, x, high, n, false);
    append_profile(&pts, x, low, n, true);
    html_open(w, "polygon");
    html_attr(w, "class", "%s", band_class[band]);
    html_attr(w, "points", "%.*s", (int)pts.len, pts.s);
    html_tag_text(w, "title", band_label[band]);
    html_close(w);
    for (size_t i = 0; i < n; i++)
      low[i] = high[i];
  }
  strbuf_free(&pts);
  free(x);
  free(low);
  free(high);
}

/* mark, under the axis, the stretches of the reference where the total
   depth <total> is below <min_depth> */
static void low_depth_marks(html_writer_t *w, const size_t *total,
    size_t len, size_t min_depth) {
  double inner = COV_W - COV_LEFT - COV_RIGHT;
  size_t from = 0;
  bool open = false;
  for (size_t p = 0; p <= len; p++) {
    bool low = p < len && total[p] < min_depth;
    if (low && !open)
      from = p;
    if (!low && open) {
      double x = COV_LEFT + (double)from * inner / (double)(len - 1);
      double width = (double)(p - from) * inner / (double)(len - 1);
      html_open(w, "rect");
      html_attr(w, "class", "c-low");
      html_attr(w, "x", "%.1f", x);
      html_attr(w, "y", "%d", COV_TOP + COV_PLOT_H + 1);
      html_attr(w, "width", "%.1f", width < 1 ? 1 : width);
      html_attr(w, "height", "3");
      html_open(w, "title");
      html_inline(w);
      html_textf(w, "%zu-%zu: below %zu\xc3\x97", from + 1, p, min_depth);
      html_close(w);
      html_close(w);
    }
    open = low;
  }
}

void plots_coverage(html_writer_t *w, const counter_t *c, size_t ref,
    size_t len, size_t min_depth, bool mark_site) {
  assert(c != NULL);
  if (len < 2)
    return;
  /* the depth by band and position, the total depth, and that of a
     class of reads */
  size_t *depth, *total, *one;
  XCALLOC(depth, N_BANDS * len, sizeof(size_t));
  XCALLOC(total, len, sizeof(size_t));
  XMALLOC(one, sizeof(size_t) * len);
  size_t max = 0;
  for (size_t k = 0; k < COVER_N_CATEGORIES; k++) {
    counter_coverage(c, ref, k, one);
    for (size_t p = 0; p < len; p++) {
      depth[(size_t)band_of_cover[k] * len + p] += one[p];
      total[p] += one[p];
    }
  }
  free(one);
  for (size_t p = 0; p < len; p++)
    if (total[p] > max)
      max = total[p];
  if (max == 0) {
    free(depth);
    free(total);
    markup_note(w, "No located reads.");
    return;
  }
  double inner = COV_W - COV_LEFT - COV_RIGHT, bottom = COV_TOP + COV_PLOT_H;
  html_open(w, "svg");
  html_attr(w, "viewBox", "0 0 %d %d", COV_W, COV_H);
  html_attr(w, "role", "img");
  coverage_bands(w, depth, len, max);
  if (min_depth > 0)
    low_depth_marks(w, total, len, min_depth);
  free(depth);
  free(total);
  svg_line(w, "axis", COV_LEFT, bottom, COV_LEFT + inner, bottom);
  svg_label(w, COV_LEFT - 4, COV_TOP + 4, "end", "%s\xc3\x97",
            html_num(max).s);
  svg_label(w, COV_LEFT, bottom + 13, "start", "1");
  svg_label(w, COV_LEFT + inner, bottom + 13, "end", "%s", html_num(len).s);
  /* the depth from which a position is well covered (not drawn when it
     would lie on the axis or on the top of the plot) */
  if (min_depth > 0 && min_depth * 20 >= max && min_depth * 10 <= max * 9) {
    double y = bottom - COV_PLOT_H * (double)min_depth / (double)max;
    svg_line(w, "mark", COV_LEFT, y, COV_LEFT + inner, y);
    svg_label(w, COV_LEFT - 4, y + 4, "end", "%zu\xc3\x97", min_depth);
  }
  size_t site;
  if (mark_site && main_site(c, ref, len, &site)) {
    double x = COV_LEFT + (double)site * inner / (double)(len - 1);
    svg_line(w, "site", x, COV_TOP - 2, x, bottom);
    bool right_half = site * 2 > len;
    html_open(w, "text");
    html_attr(w, "class", "sitel");
    html_attr(w, "x", "%.1f", right_half ? x - 3 : x + 3);
    html_attr(w, "y", "%d", COV_TOP + 6);
    html_attr(w, "text-anchor", "%s", right_half ? "end" : "start");
    html_inline(w);
    html_textf(w, "%zu", site + 1);
    html_close(w);
  }
  html_close(w);
}

/* --- lengths of the complete amplicons --- */

/* geometry of the plot of the lengths (pixels of the viewBox) */
enum { LEN_W = 300, LEN_H = 122, LEN_LEFT = 44, LEN_RIGHT = 8, LEN_TOP = 12,
       LEN_PLOT_H = 92, LEN_BARS = 66 };

/* the end of the x axis: a little beyond the reference, or, without one,
   the length below which nearly all the amplicons are */
static size_t lengths_axis_end(const lengths_t *l, size_t group,
    size_t ref_len) {
  if (ref_len > 0)
    return ref_len + ref_len / 12;
  size_t total = 0, sum = 0, end = LENGTHS_BIN;
  for (size_t b = 0; b < LENGTHS_N_BINS; b++)
    total += lengths_count(l, group, b);
  for (size_t b = 0; b < LENGTHS_N_BINS && sum * 1000 < total * 995; b++) {
    sum += lengths_count(l, group, b);
    end = (b + 1) * LENGTHS_BIN;
  }
  return end;
}

void plots_lengths(html_writer_t *w, const lengths_t *l, size_t group,
    size_t ref_len) {
  assert(l != NULL);
  size_t end = lengths_axis_end(l, group, ref_len);
  /* the bins of the lengths, gathered in the bars (the amplicons longer
     than the axis in the last one) */
  size_t bars[LEN_BARS] = {0}, max = 0;
  for (size_t b = 0; b < LENGTHS_N_BINS; b++) {
    size_t i = b * LENGTHS_BIN * LEN_BARS / end;
    bars[i < LEN_BARS ? i : LEN_BARS - 1] += lengths_count(l, group, b);
  }
  for (size_t i = 0; i < LEN_BARS; i++)
    if (bars[i] > max)
      max = bars[i];
  if (max == 0) {
    markup_note(w, "No complete amplicons.");
    return;
  }
  double inner = LEN_W - LEN_LEFT - LEN_RIGHT, bottom = LEN_TOP + LEN_PLOT_H;
  double bar_w = inner / LEN_BARS;
  html_open(w, "svg");
  html_attr(w, "viewBox", "0 0 %d %d", LEN_W, LEN_H);
  html_attr(w, "role", "img");
  for (size_t i = 0; i < LEN_BARS; i++) {
    if (bars[i] == 0)
      continue;
    /* (a bar with amplicons is at least one pixel high) */
    double h = LEN_PLOT_H * (double)bars[i] / (double)max;
    size_t from = i * end / LEN_BARS, to = (i + 1) * end / LEN_BARS;
    bool is_short = ref_len > 0 &&
                    (double)to <= LENGTHS_SHORT_FRACTION * (double)ref_len;
    html_open(w, "rect");
    html_attr(w, "class", "%s", is_short ? "c-short" : "c-full");
    html_attr(w, "x", "%.1f", LEN_LEFT + (double)i * bar_w);
    html_attr(w, "y", "%.1f", bottom - (h < 1 ? 1 : h));
    html_attr(w, "width", "%.1f", bar_w - 0.6);
    html_attr(w, "height", "%.1f", h < 1 ? 1 : h);
    html_open(w, "title");
    html_inline(w);
    html_textf(w, "%zu-%zu bases: %s", from, to - 1, html_num(bars[i]).s);
    html_close(w);
    html_close(w);
  }
  svg_line(w, "axis", LEN_LEFT, bottom, LEN_LEFT + inner, bottom);
  svg_label(w, LEN_LEFT - 4, LEN_TOP + 4, "end", "%s", html_num(max).s);
  svg_label(w, LEN_LEFT, bottom + 13, "start", "0");
  if (ref_len > 0) {
    double x = LEN_LEFT + inner * (double)ref_len / (double)end;
    svg_line(w, "mark", x, LEN_TOP - 2, x, bottom);
    svg_label(w, x, bottom + 13, "middle", "%s", html_num(ref_len).s);
  } else
    svg_label(w, LEN_LEFT + inner, bottom + 13, "end", "%s",
              html_num(end).s);
  html_close(w);
}

/* --- median depth of the segments in the samples --- */

/* geometry of the plot of the depths (pixels of the viewBox) */
enum { DEP_W = 1000, DEP_H = 190, DEP_LEFT = 52, DEP_RIGHT = 36,
       DEP_TOP = 10, DEP_PLOT_H = 150 };

/* the height of the depth <depth> above the axis, on a logarithmic scale
   up to 10^<decades> (the depths below 1 are on the axis) */
static double depth_height(double depth, int decades) {
  return depth < 1 ? 0 : DEP_PLOT_H * log10(depth) / decades;
}

void plots_segment_depths(html_writer_t *w, const size_t *depths,
    const char *const *names, const bool *judged, size_t n,
    size_t min_depth) {
  assert(depths != NULL);
  assert(names != NULL);
  size_t max = 1000;
  for (size_t i = 0; i < n * SEGMENTS_N; i++)
    if (depths[i] > max)
      max = depths[i];
  int decades = (int)ceil(log10((double)max));
  double inner = DEP_W - DEP_LEFT - DEP_RIGHT, bottom = DEP_TOP + DEP_PLOT_H;
  double column = inner / SEGMENTS_N;
  html_open(w, "svg");
  html_attr(w, "viewBox", "0 0 %d %d", DEP_W, DEP_H);
  html_attr(w, "role", "img");
  /* (a band behind every other segment keeps their dots apart) */
  for (int s = 1; s <= SEGMENTS_N; s += 2) {
    html_elem(w, "rect");
    html_attr(w, "class", "band");
    html_attr(w, "x", "%.1f", DEP_LEFT + column * (s - 1));
    html_attr(w, "y", "%d", DEP_TOP);
    html_attr(w, "width", "%.1f", column);
    html_attr(w, "height", "%d", DEP_PLOT_H);
  }
  double tick = 1;
  for (int d = 0; d <= decades; d++, tick *= 10) {
    double y = bottom - depth_height(tick, decades);
    svg_line(w, d == 0 ? "axis" : "grid", DEP_LEFT, y, DEP_LEFT + inner, y);
    svg_label(w, DEP_LEFT - 5, y + 4, "end", "%s%s",
              d == 0 ? "\xe2\x89\xa4" : "", html_num((size_t)tick).s);
  }
  if (min_depth > 0) {
    double y = bottom - depth_height((double)min_depth, decades);
    svg_line(w, "mark", DEP_LEFT, y, DEP_LEFT + inner, y);
    svg_label(w, DEP_LEFT + inner + 4, y + 4, "start", "%zu\xc3\x97",
              min_depth);
  }
  for (int s = 1; s <= SEGMENTS_N; s++) {
    double center = DEP_LEFT + column * (s - 0.5);
    svg_label(w, center, bottom + 16, "middle", "%s", segments_gene(s));
    for (size_t i = 0; i < n; i++) {
      size_t depth = depths[i * SEGMENTS_N + (size_t)s - 1];
      /* (the samples side by side, in their order, within the column) */
      double x = center + (n > 1 ? column * 0.7 *
                 ((double)i / (double)(n - 1) - 0.5) : 0);
      html_open(w, "circle");
      html_attr(w, "class", judged[i] ? "dot" : "dot open");
      html_attr(w, "cx", "%.1f", x);
      html_attr(w, "cy", "%.1f",
                bottom - depth_height((double)depth, decades));
      html_attr(w, "r", "%s", n > 24 ? "2.5" : "4");
      html_open(w, "title");
      html_inline(w);
      html_textf(w, "%s, %s: %s\xc3\x97", names[i], segments_gene(s),
                 html_num(depth).s);
      html_close(w);
      html_close(w);
    }
  }
  html_close(w);
}
