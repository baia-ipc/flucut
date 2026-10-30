/* (c) 2026 Giorgio Gonnella */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "file.h"
#include "lengths.h"
#include "logger.h"
#include "results.h"
#include "segments.h"
#include "strbuf.h"

/* the peaks of the lengths shown without references */
#define MAX_PEAKS 12

/* the header and a line of the summary table */
static void summary_header(FILE *f) {
  fputs("output\treads_in\tbases_in\treads_out\tbases_out\t"
        "both_ends\tone_end\tno_end\tminus_strand\tdimers\ttoo_short\t"
        "with_internal_match\tleft_terminal\tleft_internal\t"
        "left_unknown\tright_terminal\tright_internal\tright_unknown\t"
        "short_complete\n", f);
}

static void summary_line(FILE *f, const char *name,
    const counter_reads_t *s, bool refs, size_t short_complete) {
  fprintf(f, "%s\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu",
          name, s->reads_in, s->bases_in, s->reads_out, s->bases_out,
          s->n_ends[2], s->n_ends[1], s->n_ends[0], s->minus, s->dimers,
          s->too_short, s->with_internal);
  /* (the origins are only known with references) */
  static const end_origin_t order[] = {ORIGIN_TERMINAL, ORIGIN_INTERNAL,
                                       ORIGIN_UNKNOWN};
  for (size_t e = 0; e < 2; e++)
    for (size_t o = 0; o < 3; o++) {
      if (refs)
        fprintf(f, "\t%zu", s->origin[e][order[o]]);
      else
        fputs("\tNA", f);
    }
  if (refs)
    fprintf(f, "\t%zu\n", short_complete);
  else
    fputs("\tNA\n", f);
}

/* log the counts of an output */
void results_log_counts(const char *name, const counter_reads_t *s, bool refs,
    Logger *logger) {
  lprintf(logger, LOG_SUCCESS, "%s: %zu of %zu reads written (%zu with "
          "both primers, %zu with one, %zu with none; %zu too short, of "
          "which %zu dimers)\n", name, s->reads_out, s->reads_in,
          s->n_ends[2], s->n_ends[1], s->n_ends[0], s->too_short,
          s->dimers);
  if (s->adapter_start + s->adapter_end > 0)
    lprintf(logger, LOG_INFO, "%s: adapters removed at the ends without "
            "a primer: %zu at the start of a read, %zu at the end\n", name,
            s->adapter_start, s->adapter_end);
  if (refs)
    lprintf(logger, LOG_INFO, "%s: ends started by internal priming: "
            "%zu left, %zu right; not located in the references: %zu "
            "left, %zu right\n", name, s->origin[0][ORIGIN_INTERNAL],
            s->origin[1][ORIGIN_INTERNAL], s->origin[0][ORIGIN_UNKNOWN],
            s->origin[1][ORIGIN_UNKNOWN]);
}

/* the name of a group of the lengths (see counter_lengths()) */
static const char* group_name(const fasta_t *refs, size_t group) {
  if (refs == NULL)
    return "all";
  return group < refs->n_records ? refs->records[group].name :
         "not_located";
}

/*
   Write the lengths of the complete amplicons of an output to <f>
   (lengths.tsv) and log the check (see lengths.h): with references, the
   groups with short amplicons; without, the peaks of the lengths.
   Returns the number of short amplicons of the references.
*/
static size_t write_lengths(const lengths_t *l, const fasta_t *refs,
    const char *output, const char *path, FILE *f, Logger *logger) {
  size_t n_short = 0;
  size_t n_groups = lengths_n_groups(l);
  for (size_t g = 0; g < n_groups; g++) {
    const char *name = group_name(refs, g);
    size_t len = strlen(output) + strlen(name) + 3;
    char *prefix;
    XMALLOC(prefix, len);
    snprintf(prefix, len, "%s\t%s\t", output, name);
    lengths_write(l, g, prefix, f);
    free(prefix);
    if (refs == NULL) {
      size_t peaks[MAX_PEAKS];
      size_t n = lengths_peaks(l, g, peaks, MAX_PEAKS);
      if (n == 0)
        continue;
      strbuf_t line;
      strbuf_init(&line);
      for (size_t i = 0; i < n; i++) {
        char peak[64];
        snprintf(peak, sizeof(peak), " %zu-%zu (%zu)",
                 peaks[i] * LENGTHS_BIN,
                 peaks[i] * LENGTHS_BIN + LENGTHS_BIN - 1,
                 lengths_count(l, g, peaks[i]));
        strbuf_append_str(&line, peak);
      }
      lprintf(logger, LOG_INFO, "%s: complete amplicons, lengths "
              "between the primers (peaks):%.*s\n", path, (int)line.len,
              line.s);
      strbuf_free(&line);
      continue;
    }
    if (g == refs->n_records)
      continue;
    lengths_summary_t s;
    lengths_summarize(l, g, refs->records[g].seq_len, &s);
    n_short += s.short_count;
    if (s.short_count == 0)
      continue;
    /* (shorter: below LENGTHS_SHORT_FRACTION of the reference) */
    int segment = segments_of_description(refs->records[g].description);
    char seg[32] = "";
    if (segment > 0)
      snprintf(seg, sizeof(seg), "%s segment: ", segments_gene(segment));
    lprintf(logger, LOG_WARN, "%s: %s%zu of %zu complete amplicons "
            "shorter than the closest reference (%s, %zu bases), most at "
            "%zu-%zu bases: internal priming or deletions\n", path, seg,
            s.short_count, s.total, name, refs->records[g].seq_len,
            s.short_bin * LENGTHS_BIN,
            s.short_bin * LENGTHS_BIN + LENGTHS_BIN - 1);
  }
  return n_short;
}

void results_write(const char *outdir, const results_group_t *groups,
    size_t n, const counter_t *total, const fasta_t *refs,
    const run_report_t *report, bool running, Logger *logger) {
  char *summary_path = path_join(outdir, RESULTS_SUMMARY_FILE);
  FILE *summary = xfopen(summary_path, REPORT_TMP_SFX, "w",
                         "summary table", logger);
  summary_header(summary);
  char *lengths_path = path_join(outdir, RESULTS_LENGTHS_FILE);
  FILE *lengths_file = xfopen(lengths_path, REPORT_TMP_SFX, "w",
                              "lengths table", logger);
  fputs("output\tgroup\tlength\treads\n", lengths_file);
  /* (in a run in progress, the messages are given at the end only) */
  Logger *log = running ? NULL : logger;
  for (size_t i = 0; i < n; i++) {
    const results_group_t *g = groups + i;
    const counter_reads_t *stats = counter_reads(g->counts);
    if (!running)
      results_log_counts(g->output, stats, refs != NULL, logger);
    size_t n_short = write_lengths(counter_lengths(g->counts), refs,
        g->name, g->output, lengths_file, log);
    summary_line(summary, g->name, stats, refs != NULL, n_short);
  }
  xfclose(&summary, logger);
  xfclose(&lengths_file, logger);
  xfrename(summary_path, REPORT_TMP_SFX, summary_path, NULL, logger);
  xfrename(lengths_path, REPORT_TMP_SFX, lengths_path, NULL, logger);
  report_section_t *sections;
  XMALLOC(sections, sizeof(report_section_t) * (n + 1));
  for (size_t i = 0; i < n; i++)
    sections[i] = (report_section_t){groups[i].name, groups[i].counts,
                                     groups[i].n_files};
  if (!running)
    lprintf(logger, LOG_PHASE, "Writing the reports...\n");
  run_report_write(report, sections, n, total, running, logger);
  free(sections);
  if (!running) {
    lprintf(logger, LOG_INFO, "Tables written: %s, %s\n", summary_path,
            lengths_path);
    lprintf(logger, LOG_SUCCESS, "Reports written: %s/%s, %s/%s\n",
            outdir, RUN_REPORT_HTML_FILE, outdir, RUN_REPORT_JSON_FILE);
  }
  free(summary_path);
  free(lengths_path);
}
