/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include "html-page.h"
#include "html-report-markup.h"
#include "html-report-reads.h"
#include "html-report-sample.h"
#include "svg-writer.h"
#include "numutils.h"

/* a row of the table: a name and the texts of the two sets of reads */
static void row(html_writer_t *w, const char *name, const char *input,
    const char *written) {
  markup_open_row(w);
  markup_td(w, "%s", name);
  markup_td(w, "%s", input);
  markup_td(w, "%s", written);
  html_close(w);
}

/* a row with a number of each set */
static void number_row(html_writer_t *w, const char *name, size_t input,
    size_t written) {
  row(w, name, html_num(input).s, html_num(written).s);
}

/* the table of the two sets of reads */
static void write_table(html_writer_t *w, const read_stats_summary_t *in,
    const read_stats_summary_t *out) {
  char a[48], b[48];
  html_open(w, "table");
  markup_header_row(w, "", "input", "written", NULL);
  number_row(w, "reads", in->reads, out->reads);
  row(w, "bases", sample_bases(in->bases).s, sample_bases(out->bases).s);
  number_row(w, "median length", in->median_len, out->median_len);
  number_row(w, "mean length", in->mean_len, out->mean_len);
  number_row(w, "N50", in->n50, out->n50);
  snprintf(a, sizeof(a), "%s\xe2\x80\x93%s", html_num(in->min_len).s,
           html_num(in->max_len).s);
  snprintf(b, sizeof(b), "%s\xe2\x80\x93%s", html_num(out->min_len).s,
           html_num(out->max_len).s);
  row(w, "shortest\xe2\x80\x93longest", a, b);
  if (in->has_quality) {
    snprintf(a, sizeof(a), "Q%.0f", in->median_quality - 0.5);
    snprintf(b, sizeof(b), "Q%.0f", out->median_quality - 0.5);
    row(w, "median quality", a, out->has_quality ? b : "\xe2\x80\x93");
    row(w, "reads \xe2\x89\xa5Q10", fraction_percent_text(in->q10).s,
        fraction_percent_text(out->q10).s);
    row(w, "reads \xe2\x89\xa5Q20", fraction_percent_text(in->q20).s,
        fraction_percent_text(out->q20).s);
  }
  html_close(w);
}

/* geometry of the plots of the distributions (pixels of the viewBox) */
enum { DIST_W = 360, DIST_H = 150, DIST_LEFT = 44, DIST_RIGHT = 10,
       DIST_TOP = 10, DIST_PLOT_H = 112, DIST_BARS = 64 };

/* The distribution of a value in the two sets of reads: the counts <in>
   and <out> of <n> bins of width <width>, up to the bin <end>
   (excluded): the input as bars, the reads written as an outline in
   steps over the same bins */
static void distribution(html_writer_t *w, const size_t *in,
    const size_t *out, size_t width, size_t end) {
  /* the bins gathered in at most DIST_BARS bars */
  size_t per_bar = (end + DIST_BARS - 1) / DIST_BARS;
  size_t n_bars = (end + per_bar - 1) / per_bar, max = 0;
  size_t bars[2][DIST_BARS] = {{0}};
  for (size_t b = 0; b < end; b++) {
    bars[0][b / per_bar] += in[b];
    bars[1][b / per_bar] += out[b];
  }
  for (size_t i = 0; i < n_bars; i++) {
    if (bars[0][i] > max)
      max = bars[0][i];
    if (bars[1][i] > max)
      max = bars[1][i];
  }
  if (max == 0)
    return;
  double inner = DIST_W - DIST_LEFT - DIST_RIGHT;
  double bottom = DIST_TOP + DIST_PLOT_H, bar_w = inner / (double)n_bars;
  html_open(w, "svg");
  html_attr(w, "viewBox", "0 0 %d %d", DIST_W, DIST_H);
  html_attr(w, "role", "img");
  for (size_t i = 0; i < n_bars; i++) {
    if (bars[0][i] == 0)
      continue;
    double h = DIST_PLOT_H * (double)bars[0][i] / (double)max;
    html_elem(w, "rect");
    html_attr(w, "class", "c-oth");
    html_attr(w, "x", "%.1f", DIST_LEFT + (double)i * bar_w);
    html_attr(w, "y", "%.1f", bottom - h);
    html_attr(w, "width", "%.1f", bar_w > 1.5 ? bar_w - 0.6 : bar_w);
    html_attr(w, "height", "%.1f", h);
  }
  /* the reads written: an outline in steps, on the edges of the bars */
  char points[DIST_BARS * 32];
  size_t len = 0;
  for (size_t i = 0; i < n_bars; i++) {
    double y = bottom - DIST_PLOT_H * (double)bars[1][i] / (double)max;
    int n = snprintf(points + len, sizeof(points) - len,
                     "%s%.1f,%.1f %.1f,%.1f", i > 0 ? " " : "",
                     DIST_LEFT + (double)i * bar_w, y,
                     DIST_LEFT + (double)(i + 1) * bar_w, y);
    if (n > 0 && len + (size_t)n < sizeof(points))
      len += (size_t)n;
  }
  html_elem(w, "polyline");
  html_attr(w, "class", "line");
  html_attr(w, "points", "%s", points);
  svg_line(w, "axis", DIST_LEFT, bottom, DIST_LEFT + inner, bottom);
  svg_label(w, DIST_LEFT - 4, DIST_TOP + 4, "end", "%s", html_num(max).s);
  svg_label(w, DIST_LEFT, bottom + 14, "start", "0");
  svg_label(w, DIST_LEFT + inner / 2, bottom + 14, "middle", "%s",
            html_num(n_bars * per_bar * width / 2).s);
  svg_label(w, DIST_LEFT + inner, bottom + 14, "end", "%s",
            html_num(n_bars * per_bar * width).s);
  html_close(w);
}

/* the bin below which nearly all (99.5%) of the <total> counts of the
   <n> bins <counts> are (at least 1) */
static size_t distribution_end(const size_t *counts, size_t n,
    size_t total) {
  size_t sum = 0, end = 1;
  for (size_t b = 0; b < n && sum * 1000 < total * 995; b++) {
    sum += counts[b];
    end = b + 1;
  }
  return end;
}

/* a plot with its title and the legend of the two sets */
static void write_plot(html_writer_t *w, const char *title,
    const size_t *in, const size_t *out, size_t n, size_t width,
    size_t total) {
  html_open(w, "div");
  html_tag_text(w, "b", title);
  markup_open(w, "div", "lg");
  markup_swatch(w, "c-oth", "input");
  markup_swatch(w, "c-line", "written");
  html_close(w);
  size_t total_out = 0;
  for (size_t b = 0; b < n; b++)
    total_out += out[b];
  size_t end = distribution_end(in, n, total);
  size_t end_out = distribution_end(out, n, total_out);
  distribution(w, in, out, width, end > end_out ? end : end_out);
  html_close(w);
}

/* what the titles of the reads say about the run */
static void write_metadata(html_writer_t *w, const read_metadata_t *m) {
  if (!read_metadata_known(m))
    return;
  markup_open(w, "dl", "kv meta");
  for (size_t i = 0; i < READ_METADATA_N; i++) {
    const char *text = read_metadata_text(m, i);
    if (text == NULL)
      continue;
    html_open(w, "div");
    html_tag_text(w, "dt", read_metadata_label(i));
    html_tag_text(w, "dd", text);
    html_close(w);
  }
  html_close(w);
}

void reads_write(html_writer_t *w, const counter_t *c) {
  assert(c != NULL);
  const read_stats_t *stats = counter_read_stats(c);
  const read_stats_set_t *in = stats->sets + READ_STATS_INPUT;
  const read_stats_set_t *out = stats->sets + READ_STATS_WRITTEN;
  read_stats_summary_t summary_in, summary_out;
  read_stats_summarize(stats, READ_STATS_INPUT, &summary_in);
  read_stats_summarize(stats, READ_STATS_WRITTEN, &summary_out);
  markup_open(w, "section", "card");
  html_tag_text(w, "h3", "Read length and quality");
  markup_open(w, "div", "readqc");
  html_open(w, "div");
  write_table(w, &summary_in, &summary_out);
  html_close(w);
  write_plot(w, "Reads by length (bases)", in->len_reads, out->len_reads,
             READ_STATS_LEN_BINS, READ_STATS_LEN_BIN, in->reads);
  if (summary_in.has_quality)
    write_plot(w, "Reads by quality (Q)", in->quality, out->quality,
               READ_STATS_QUAL_BINS, 1, in->with_quality);
  html_close(w);
  write_metadata(w, counter_metadata(c));
  html_close(w);
}
