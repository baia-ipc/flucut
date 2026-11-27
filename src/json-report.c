/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "file.h"
#include "json-coverage.h"
#include "json-reads.h"
#include "json-report.h"
#include "json-sample.h"
#include "json-writer.h"
#include "qc-summary.h"
#include "segments.h"
#include "virus.h"
#include "strutils.h"

static void write_reads(json_writer_t *w, const counter_reads_t *s) {
  json_member(w, "reads");
  json_open_block(w, '{');
  json_count_member(w, "in", s->reads_in);
  json_count_member(w, "bases_in", s->bases_in);
  json_count_member(w, "written", s->reads_out);
  json_count_member(w, "bases_written", s->bases_out);
  json_count_member(w, "both_primers", s->n_ends[2]);
  json_count_member(w, "one_primer", s->n_ends[1]);
  json_count_member(w, "no_primer", s->n_ends[0]);
  json_count_member(w, "dimers", s->dimers);
  json_count_member(w, "too_short", s->too_short);
  json_count_member(w, "minus_strand", s->minus);
  json_count_member(w, "with_internal_match", s->with_internal);
  json_count_member(w, "split", s->split);
  json_count_member(w, "cut_at_junction", s->cut);
  json_count_member(w, "adapter_start", s->adapter_start);
  json_count_member(w, "adapter_end", s->adapter_end);
  json_close_block(w, '}', false);
}

/* the reads written, by category (null without references) */
static void write_categories(json_writer_t *w, const report_t *r,
    const counter_reads_t *s) {
  json_member(w, "written_by_category");
  if (r->references == NULL) {
    fputs("null", w->f);
    return;
  }
  fputc('{', w->f);
  bool first = true;
  for (size_t k = 0; k < READ_N_CATEGORIES; k++) {
    json_sep(w->f, &first);
    fprintf(w->f, "\"%s\": %zu", read_category_name(k), s->written[k]);
  }
  fputc('}', w->f);
}

/* the primers found at each end, by primer of the scheme */
static void write_primers(json_writer_t *w, const report_t *r,
    const counter_t *c) {
  json_member(w, "primers");
  json_open_block(w, '{');
  for (size_t e = 0; e < COUNTER_N_ENDS; e++) {
    const primer_t *p = e == COUNTER_LEFT ? r->scheme->left :
                        r->scheme->right;
    json_member(w, report_end_name(e));
    fputc('{', w->f);
    bool first = true;
    for (size_t v = 0; v < counter_n_variants(c, e); v++) {
      json_sep(w->f, &first);
      fputs_json_string(p[v].name, w->f);
      fprintf(w->f, ": %zu", counter_primer(c, e, v));
    }
    fputc('}', w->f);
  }
  json_close_block(w, '}', false);
}

/* the read ends by origin and the internal priming sites (null without
   references) */
static void write_internal(json_writer_t *w, const report_t *r,
    const counter_t *c) {
  FILE *f = w->f;
  json_member(w, "ends");
  if (r->references == NULL)
    fputs("null", f);
  else {
    const counter_reads_t *s = counter_reads(c);
    json_open_block(w, '{');
    for (size_t e = 0; e < COUNTER_N_ENDS; e++) {
      json_member(w, report_end_name(e));
      fprintf(f, "{\"terminal\": %zu, \"internal\": %zu, \"unknown\": "
              "%zu}", s->origin[e][ORIGIN_TERMINAL],
              s->origin[e][ORIGIN_INTERNAL], s->origin[e][ORIGIN_UNKNOWN]);
    }
    json_close_block(w, '}', false);
  }
  json_member(w, "internal_sites");
  if (r->references == NULL) {
    fputs("null", f);
    return;
  }
  size_t n;
  report_site_t *sites = report_sites(r->references, c, &n);
  json_open_block(w, '[');
  for (size_t i = 0; i < n && i < JSON_MAX_SITES; i++) {
    const report_site_t *x = sites + i;
    json_member(w, NULL);
    fputs("{\"reference\": ", f);
    fputs_json_string(r->references->records[x->ref].name, f);
    fprintf(f, ", \"end\": \"%s\", \"position\": %zu, \"from\": %zu, "
            "\"to\": %zu, \"reads\": %zu}", report_end_name(x->end),
            x->pos, x->from, x->to, x->reads);
  }
  json_close_block(w, ']', true);
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
  json_member(w, "short_complete");
  if (r->references == NULL)
    fputs("null", f);
  else
    fprintf(f, "%zu", report_short_complete(r->references, c));
  json_member(w, "complete_amplicons");
  json_open_block(w, '[');
  for (size_t g = 0; g < lengths_n_groups(l); g++) {
    report_lengths_t x;
    report_lengths(r->references, c, g, &x);
    if (x.summary.total == 0)
      continue;
    bool ref = x.ref_len > 0;
    json_member(w, NULL);
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
      json_sep(f, &first);
      fprintf(f, "%zu", lengths_count(l, g, b));
    }
    fputs("]}", f);
  }
  json_close_block(w, ']', true);
}

/*
   What the reads say about the cleavage site of the haemagglutinin
   (cleavage-site.h); null when it was not read (no HA reference with
   the fusion peptide, or no reference at all)
*/
static void write_cleavage_site(json_writer_t *w, const report_t *r,
    const counter_t *c) {
  FILE *f = w->f;
  json_member(w, "cleavage_site");
  cleavage_result_t s;
  report_cleavage(r, c, &s);
  if (s.searched == 0) {
    fputs("null", f);
    return;
  }
  json_open_block(w, '{');
  json_member(w, "chemistry");
  fprintf(f, "\"%s\"", cleavage_chemistry_name(s.chemistry));
  json_count_member(w, "reads_searched", s.searched);
  json_count_member(w, "reads", s.reads);
  json_count_member(w, "distinct_motifs", s.distinct);
  json_count_member(w, "distinct_sites", s.distinct_proteins);
  json_count_member(w, "supported_sites", s.supported_proteins);
  json_member(w, "motif");
  if (!s.found) {
    fputs("null", f);
    json_member(w, "protein");
    fputs("null", f);
  } else {
    fputs_json_string(s.nt, f);
    json_member(w, "protein");
    fputs_json_string(s.aa, f);
  }
  json_count_member(w, "motif_reads", s.motif_reads);
  json_share_member(w, "support", s.found ? s.support : 0.0);
  json_share_member(w, "multibasic_share", s.multibasic_share);
  /* what the same reads agree on beside the site, where nothing is
     difficult: reported, not used to judge the site */
  json_count_member(w, "control_reads", s.control_reads);
  json_share_member(w, "control_support", s.control_support);
  json_member(w, "reads_disagree_beside_the_site");
  fputs(s.control_weak ? "true" : "false", f);
  json_count_member(w, "basic_residues", s.basic);
  json_member(w, "reads_as");
  fprintf(f, "\"%s\"", cleavage_reading_name(s.reading));
  /* the fusion peptide: the anchor, and what the reads carry there */
  json_member(w, "fusion_peptide");
  if (s.fusion_aa[0] == '\0')
    fputs("null", f);
  else
    fputs_json_string(s.fusion_aa, f);
  json_count_member(w, "fusion_peptide_reads", s.fusion_reads);
  json_share_member(w, "fusion_peptide_support", s.fusion_support);
  json_member(w, "reference_fusion_peptide");
  if (s.reference_fusion[0] == '\0')
    fputs("null", f);
  else
    fputs_json_string(s.reference_fusion, f);
  json_member(w, "reference_site");
  if (s.reference_aa[0] == '\0')
    fputs("null", f);
  else
    fputs_json_string(s.reference_aa, f);
  json_member(w, "disagrees_with_reference");
  fputs(s.disagrees_with_reference ? "true" : "false", f);
  json_member(w, "verdict");
  fprintf(f, "\"%s\"", cleavage_verdict_name(s.verdict));
  /* the sites the reads carry, most frequent first: more than one means
     the reads of more than one haemagglutinin are in the sample */
  json_member(w, "sites");
  json_open_block(w, '[');
  for (size_t i = 0; i < s.n_top; i++) {
    json_member(w, NULL);
    fputs("{\"protein\": ", f);
    fputs_json_string(s.top[i].aa, f);
    fputs(", \"motif\": ", f);
    fputs_json_string(s.top[i].nt, f);
    fprintf(f, ", \"reads\": %zu, \"share\": %.4f, \"reads_as\": \"%s\"}",
            s.top[i].reads, (double)s.top[i].reads / (double)s.reads,
            cleavage_reading_name(s.top[i].reading));
  }
  json_close_block(w, ']', true);
  json_close_block(w, '}', true);
}

/* the statistics of a counter, as members of the open object */
static void write_section(json_writer_t *w, const report_t *r,
    const counter_t *c) {
  write_reads(w, counter_reads(c));
  write_categories(w, r, counter_reads(c));
  json_member(w, "read_lengths");
  json_reads_lengths(w->f, c);
  json_member(w, "read_quality");
  json_reads_quality(w->f, c);
  json_member(w, "sequencing");
  json_reads_sequencing(w->f, c);
  json_sample_write(w, r, c);
  write_cleavage_site(w, r, c);
  write_primers(w, r, c);
  write_internal(w, r, c);
  write_lengths(w, r, c);
  json_member(w, "coverage");
  if (r->references == NULL) {
    fputs("null", w->f);
    return;
  }
  json_open_block(w, '[');
  for (size_t g = 0; g < r->references->n_records; g++)
    if (json_coverage_has(c, g, r->references)) {
      json_member(w, NULL);
      json_coverage_write(w->f, r->references, c, g);
    }
  json_close_block(w, ']', true);
}

/* a member whose value is a segment (1 to 8), or null if 0 (all) */
static void segment_member(json_writer_t *w, const char *name,
    int segment) {
  json_member(w, name);
  if (segment == 0)
    fputs("null", w->f);
  else
    fprintf(w->f, "%d", segment);
}

/* the parameters of the run */
static void write_parameters(json_writer_t *w, const report_t *r) {
  FILE *f = w->f;
  const options_t *o = r->options;
  json_member(w, "scheme");
  json_open_block(w, '{');
  json_member(w, "name");
  fputs_json_string(r->scheme->name, f);
  for (size_t e = 0; e < COUNTER_N_ENDS; e++) {
    const primer_t *p = e == COUNTER_LEFT ? r->scheme->left :
                        r->scheme->right;
    size_t n = e == COUNTER_LEFT ? r->scheme->n_left : r->scheme->n_right;
    json_member(w, report_end_name(e));
    fputc('[', f);
    bool first = true;
    for (size_t v = 0; v < n; v++) {
      json_sep(f, &first);
      fputs_json_string(p[v].name, f);
    }
    fputc(']', f);
  }
  json_close_block(w, '}', false);
  json_member(w, "trim_policy");
  fputs(o->policy == TRIM_PRIMER ? "\"primer\"" : "\"terminus\"", f);
  json_member(w, "references");
  if (r->references == NULL)
    fputs("null", f);
  else {
    json_open_block(w, '{');
    json_member(w, "file");
    fputs_json_string(o->references, f);
    json_member(w, "sequences");
    json_open_block(w, '[');
    for (size_t i = 0; i < r->references->n_records; i++) {
      json_member(w, NULL);
      fputs("{\"name\": ", f);
      fputs_json_string(r->references->records[i].name, f);
      const char *description = r->references->records[i].description;
      int segment = segments_of_description(description);
      fprintf(f, ", \"length\": %zu, \"virus\": \"%s\", \"segment\": ",
              r->references->records[i].seq_len,
              virus_type_name(virus_type_of_description(description)));
      if (segment == 0)
        fputs("null, \"gene\": null}", f);
      else
        /* (the number of the segment in the type of the reference) */
        fprintf(f, "%d, \"gene\": \"%s\"}",
                segments_number(segment,
                    virus_type_of_description(description)),
                segments_gene(segment));
    }
    json_close_block(w, ']', false);
    json_close_block(w, '}', false);
  }
  json_member(w, "settings");
  json_open_block(w, '{');
  json_member(w, "max_errors");
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
  json_member(w, "platform");
  fputs(o->illumina ? "\"illumina\"" : "\"nanopore\"", f);
  json_count_member(w, "window", (size_t)o->window);
  json_count_member(w, "min_length", o->min_length);
  json_count_member(w, "probe_length", o->probe_len);
  json_member(w, "probe_errors");
  fprintf(f, "%g", o->probe_error_rate);
  json_count_member(w, "terminal_window", (size_t)o->terminal_window);
  json_member(w, "exhaustive_search");
  fputs(o->exhaustive_search ? "true" : "false", f);
  segment_member(w, "segment", o->segment);
  segment_member(w, "report_segment", o->report_segment);
  json_member(w, "tags");
  fputs(o->tags ? "true" : "false", f);
  json_member(w, "split");
  fputs(o->split ? "true" : "false", f);
  json_member(w, "by_category");
  fputs(o->by_category ? "true" : "false", f);
  json_member(w, "cut_deletions");
  fputs(o->cut_deletions ? "true" : "false", f);
  json_count_member(w, "threads", o->threads);
  json_count_member(w, "compress_level", (size_t)o->compress_level);
  json_close_block(w, '}', false);
}

static void write_summary(json_writer_t *w, const report_t *r) {
  FILE *f = w->f;
  json_open_block(w, '{');
  json_member(w, "format");
  fprintf(f, "%d", JSON_SUMMARY_FORMAT);
  json_member(w, "program");
  fputs("\"flucut\"", f);
  json_member(w, "version");
  fputs_json_string(r->version, f);
  time_t now = time(NULL);
  struct tm tm;
  char date[64] = "";
  if (localtime_r(&now, &tm) == NULL ||
      strftime(date, sizeof(date), "%Y-%m-%dT%H:%M:%S%z", &tm) == 0)
    date[0] = '\0';
  json_member(w, "updated");
  fputs_json_string(date, f);
  json_member(w, "running");
  fputs(r->running ? "true" : "false", f);
  json_member(w, "command");
  fputs_json_string(r->command_line, f);
  write_parameters(w, r);
  json_member(w, "length_bin");
  fprintf(f, "%d", LENGTHS_BIN);
  json_member(w, "length_max");
  fprintf(f, "%d", LENGTHS_MAX);
  json_member(w, "short_fraction");
  fprintf(f, "%g", LENGTHS_SHORT_FRACTION);
  json_member(w, "min_depth");
  fprintf(f, "%d", QC_MIN_DEPTH);
  json_member(w, "complete_fraction");
  fprintf(f, "%g", QC_COMPLETE_FRACTION);
  json_member(w, "all");
  json_open_block(w, '{');
  write_section(w, r, r->total);
  json_close_block(w, '}', false);
  json_member(w, "groups");
  json_open_block(w, '[');
  for (size_t i = 0; i < r->n_sections; i++) {
    const report_section_t *g = r->sections + i;
    json_member(w, NULL);
    json_open_block(w, '{');
    json_member(w, "name");
    fputs_json_string(g->name, f);
    json_count_member(w, "files", g->n_files);
    json_member(w, "processed");
    fputs(g->counts != NULL ? "true" : "false", f);
    if (g->counts != NULL)
      write_section(w, r, g->counts);
    json_close_block(w, '}', false);
  }
  json_close_block(w, ']', true);
  json_close_block(w, '}', false);
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
