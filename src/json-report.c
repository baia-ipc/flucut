/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "file.h"
#include "json-report.h"
#include "segments.h"
#include "strutils.h"

/*
   The summary is written as nested objects and arrays (blocks) with one
   member on each line, indented by the nesting level; the values of the
   members are written on their line (lists on one line: see sep()). The
   writer keeps, for each open block, whether it has a member yet, for
   the commas between the members (as in Lumipore).
*/
#define MAX_NESTING 8

typedef struct {
  FILE *f;
  int level;  /* number of open blocks */
  bool empty[MAX_NESTING];
} json_writer_t;

static void indent(FILE *f, int level) {
  for (int i = 0; i < level; i++)
    fputs("  ", f);
}

/* open a block: <c> is { or [ */
static void open_block(json_writer_t *w, char c) {
  assert(w->level < MAX_NESTING);
  fputc(c, w->f);
  w->empty[w->level++] = true;
}

/* start a member of the open block: on a new line (after a comma, unless
   the first); with "<k>": (k != NULL: object) */
static void member(json_writer_t *w, const char *k) {
  assert(w->level > 0);
  fputs(w->empty[w->level - 1] ? "\n" : ",\n", w->f);
  w->empty[w->level - 1] = false;
  indent(w->f, w->level);
  if (k != NULL) {
    fputs_json_string(k, w->f);
    fputs(": ", w->f);
  }
}

/* close the open block: <c> is } or ]; on a new line, or, if empty and
   <compact>, as {} or [] */
static void close_block(json_writer_t *w, char c, bool compact) {
  assert(w->level > 0);
  w->level--;
  if (!w->empty[w->level] || !compact) {
    fputc('\n', w->f);
    indent(w->f, w->level);
  }
  fputc(c, w->f);
}

/* the separator of the elements of a list on one line (none before the
   first: <first>, then set to false) */
static void sep(FILE *f, bool *first) {
  if (!*first)
    fputs(", ", f);
  *first = false;
}

/* a member with a count */
static void count_member(json_writer_t *w, const char *k, size_t n) {
  member(w, k);
  fprintf(w->f, "%zu", n);
}

static void write_reads(json_writer_t *w, const counter_reads_t *s) {
  member(w, "reads");
  open_block(w, '{');
  count_member(w, "in", s->reads_in);
  count_member(w, "bases_in", s->bases_in);
  count_member(w, "written", s->reads_out);
  count_member(w, "bases_written", s->bases_out);
  count_member(w, "both_primers", s->n_ends[2]);
  count_member(w, "one_primer", s->n_ends[1]);
  count_member(w, "no_primer", s->n_ends[0]);
  count_member(w, "dimers", s->dimers);
  count_member(w, "too_short", s->too_short);
  count_member(w, "minus_strand", s->minus);
  count_member(w, "with_internal_match", s->with_internal);
  count_member(w, "adapter_start", s->adapter_start);
  count_member(w, "adapter_end", s->adapter_end);
  close_block(w, '}', false);
}

/* the primers found at each end, by primer of the scheme */
static void write_primers(json_writer_t *w, const report_t *r,
    const counter_t *c) {
  member(w, "primers");
  open_block(w, '{');
  for (size_t e = 0; e < COUNTER_N_ENDS; e++) {
    const primer_t *p = e == COUNTER_LEFT ? r->scheme->left :
                        r->scheme->right;
    member(w, report_end_name(e));
    fputc('{', w->f);
    bool first = true;
    for (size_t v = 0; v < counter_n_variants(c, e); v++) {
      sep(w->f, &first);
      fputs_json_string(p[v].name, w->f);
      fprintf(w->f, ": %zu", counter_primer(c, e, v));
    }
    fputc('}', w->f);
  }
  close_block(w, '}', false);
}

/* the read ends by origin and the internal priming sites (null without
   references) */
static void write_internal(json_writer_t *w, const report_t *r,
    const counter_t *c) {
  FILE *f = w->f;
  member(w, "ends");
  if (r->references == NULL)
    fputs("null", f);
  else {
    const counter_reads_t *s = counter_reads(c);
    open_block(w, '{');
    for (size_t e = 0; e < COUNTER_N_ENDS; e++) {
      member(w, report_end_name(e));
      fprintf(f, "{\"terminal\": %zu, \"internal\": %zu, \"unknown\": "
              "%zu}", s->origin[e][ORIGIN_TERMINAL],
              s->origin[e][ORIGIN_INTERNAL], s->origin[e][ORIGIN_UNKNOWN]);
    }
    close_block(w, '}', false);
  }
  member(w, "internal_sites");
  if (r->references == NULL) {
    fputs("null", f);
    return;
  }
  size_t n;
  report_site_t *sites = report_sites(r->references, c, &n);
  open_block(w, '[');
  for (size_t i = 0; i < n && i < JSON_MAX_SITES; i++) {
    const report_site_t *x = sites + i;
    member(w, NULL);
    fputs("{\"reference\": ", f);
    fputs_json_string(r->references->records[x->ref].name, f);
    fprintf(f, ", \"end\": \"%s\", \"position\": %zu, \"from\": %zu, "
            "\"to\": %zu, \"reads\": %zu}", report_end_name(x->end),
            x->pos, x->from, x->to, x->reads);
  }
  close_block(w, ']', true);
  free(sites);
}

/* a bin start, or null */
static void bin_or_null(FILE *f, bool known, size_t bin) {
  if (known)
    fprintf(f, "%zu", bin * LENGTHS_BIN);
  else
    fputs("null", f);
}

/* the complete amplicons of each group of the lengths */
static void write_lengths(json_writer_t *w, const report_t *r,
    const counter_t *c) {
  FILE *f = w->f;
  const lengths_t *l = counter_lengths(c);
  member(w, "short_complete");
  if (r->references == NULL)
    fputs("null", f);
  else
    fprintf(f, "%zu", report_short_complete(r->references, c));
  member(w, "complete_amplicons");
  open_block(w, '[');
  for (size_t g = 0; g < lengths_n_groups(l); g++) {
    report_lengths_t x;
    report_lengths(r->references, c, g, &x);
    if (x.summary.total == 0)
      continue;
    bool ref = x.ref_len > 0;
    member(w, NULL);
    fputs("{\"group\": ", f);
    fputs_json_string(x.name, f);
    fputs(", \"reference_length\": ", f);
    if (ref)
      fprintf(f, "%zu", x.ref_len);
    else
      fputs("null", f);
    fprintf(f, ", \"amplicons\": %zu, \"main_length\": ", x.summary.total);
    bin_or_null(f, true, x.summary.main_bin);
    fputs(", \"short\": ", f);
    if (ref)
      fprintf(f, "%zu", x.summary.short_count);
    else
      fputs("null", f);
    fputs(", \"short_main_length\": ", f);
    bin_or_null(f, ref && x.summary.short_count > 0, x.summary.short_bin);
    /* the histogram up to the last non-empty bin */
    size_t last = 0;
    for (size_t b = 0; b < LENGTHS_N_BINS; b++)
      if (lengths_count(l, g, b) > 0)
        last = b + 1;
    fputs(", \"histogram\": [", f);
    bool first = true;
    for (size_t b = 0; b < last; b++) {
      sep(f, &first);
      fprintf(f, "%zu", lengths_count(l, g, b));
    }
    fputs("]}", f);
  }
  close_block(w, ']', true);
}

/* the statistics of a counter, as members of the open object */
static void write_section(json_writer_t *w, const report_t *r,
    const counter_t *c) {
  write_reads(w, counter_reads(c));
  write_primers(w, r, c);
  write_internal(w, r, c);
  write_lengths(w, r, c);
}

/* the parameters of the run */
static void write_parameters(json_writer_t *w, const report_t *r) {
  FILE *f = w->f;
  const options_t *o = r->options;
  member(w, "scheme");
  open_block(w, '{');
  member(w, "name");
  fputs_json_string(r->scheme->name, f);
  for (size_t e = 0; e < COUNTER_N_ENDS; e++) {
    const primer_t *p = e == COUNTER_LEFT ? r->scheme->left :
                        r->scheme->right;
    size_t n = e == COUNTER_LEFT ? r->scheme->n_left : r->scheme->n_right;
    member(w, report_end_name(e));
    fputc('[', f);
    bool first = true;
    for (size_t v = 0; v < n; v++) {
      sep(f, &first);
      fputs_json_string(p[v].name, f);
    }
    fputc(']', f);
  }
  close_block(w, '}', false);
  member(w, "trim_policy");
  fputs(o->policy == TRIM_PRIMER ? "\"primer\"" : "\"terminus\"", f);
  member(w, "references");
  if (r->references == NULL)
    fputs("null", f);
  else {
    open_block(w, '{');
    member(w, "file");
    fputs_json_string(o->references, f);
    member(w, "sequences");
    open_block(w, '[');
    for (size_t i = 0; i < r->references->n_records; i++) {
      member(w, NULL);
      fputs("{\"name\": ", f);
      fputs_json_string(r->references->records[i].name, f);
      int segment = segments_of_description(
          r->references->records[i].description);
      fprintf(f, ", \"length\": %zu, \"segment\": ",
              r->references->records[i].seq_len);
      if (segment == 0)
        fputs("null, \"gene\": null}", f);
      else
        fprintf(f, "%d, \"gene\": \"%s\"}", segment,
                segments_gene(segment));
    }
    close_block(w, ']', false);
    close_block(w, '}', false);
  }
  member(w, "settings");
  open_block(w, '{');
  member(w, "max_errors");
  const report_setting_t *s = r->max_errors;
  if (s == NULL || !s->automatic)
    fprintf(f, "{\"mode\": \"fixed\", \"value\": %g}",
            s == NULL ? o->max_errors : s->value);
  else if (s->pending)
    fputs("{\"mode\": \"auto\", \"value\": null, \"error_rate\": null, "
          "\"reads\": null, \"file\": null}", f);
  else {
    fprintf(f, "{\"mode\": \"auto\", \"value\": %g, \"error_rate\": %.4f, "
            "\"reads\": %zu, \"file\": ", s->value, s->measured,
            s->n_reads);
    fputs_json_string(s->file, f);
    fputc('}', f);
  }
  member(w, "platform");
  fputs(o->illumina ? "\"illumina\"" : "\"nanopore\"", f);
  count_member(w, "window", (size_t)o->window);
  count_member(w, "min_length", o->min_length);
  count_member(w, "probe_length", o->probe_len);
  member(w, "probe_errors");
  fprintf(f, "%g", o->probe_error_rate);
  count_member(w, "terminal_window", (size_t)o->terminal_window);
  member(w, "exhaustive_search");
  fputs(o->exhaustive_search ? "true" : "false", f);
  member(w, "tags");
  fputs(o->tags ? "true" : "false", f);
  count_member(w, "threads", o->threads);
  count_member(w, "compress_level", (size_t)o->compress_level);
  close_block(w, '}', false);
}

static void write_summary(json_writer_t *w, const report_t *r) {
  FILE *f = w->f;
  open_block(w, '{');
  member(w, "format");
  fprintf(f, "%d", JSON_SUMMARY_FORMAT);
  member(w, "program");
  fputs("\"flucut\"", f);
  member(w, "version");
  fputs_json_string(r->version, f);
  time_t now = time(NULL);
  struct tm tm;
  char date[64] = "";
  if (localtime_r(&now, &tm) == NULL ||
      strftime(date, sizeof(date), "%Y-%m-%dT%H:%M:%S%z", &tm) == 0)
    date[0] = '\0';
  member(w, "updated");
  fputs_json_string(date, f);
  member(w, "running");
  fputs(r->running ? "true" : "false", f);
  member(w, "command");
  fputs_json_string(r->command_line, f);
  write_parameters(w, r);
  member(w, "length_bin");
  fprintf(f, "%d", LENGTHS_BIN);
  member(w, "length_max");
  fprintf(f, "%d", LENGTHS_MAX);
  member(w, "short_fraction");
  fprintf(f, "%g", LENGTHS_SHORT_FRACTION);
  member(w, "all");
  open_block(w, '{');
  write_section(w, r, r->total);
  close_block(w, '}', false);
  member(w, "groups");
  open_block(w, '[');
  for (size_t i = 0; i < r->n_sections; i++) {
    const report_section_t *g = r->sections + i;
    member(w, NULL);
    open_block(w, '{');
    member(w, "name");
    fputs_json_string(g->name, f);
    count_member(w, "files", g->n_files);
    member(w, "processed");
    fputs(g->counts != NULL ? "true" : "false", f);
    if (g->counts != NULL)
      write_section(w, r, g->counts);
    close_block(w, '}', false);
  }
  close_block(w, ']', true);
  close_block(w, '}', false);
  fputc('\n', f);
}

void json_report_write(const char *path, const report_t *report,
    Logger *logger) {
  assert(path != NULL);
  assert(report != NULL);
  assert(report->total != NULL);
  FILE *f = xfopen(path, REPORT_TMP_SFX, "w", "JSON summary", logger);
  json_writer_t w = {.f = f};
  write_summary(&w, report);
  xfclose(&f, logger);
  xfrename(path, REPORT_TMP_SFX, path, NULL, logger);
}
