/* (c) 2026 Giorgio Gonnella */

/*
   flucut: removal of the primers of the influenza A multisegment RT-PCR
   from the reads, before the assembly (see README.md).

   Each group of input files (input-groups.h: a file of reads, a
   directory of files of reads, or each barcode directory of a MinKNOW
   fastq_pass directory) gives one gzipped FASTQ file in the output
   directory, and a line of the summary table (summary.tsv). The reads
   are read in batches; each batch is split between the threads, each of
   which trims and compresses its part; the parts are written in the
   order of the input, so that the output does not depend on the number
   of threads.
*/

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "banner.h"
#include "counter.h"
#include "fasta.h"
#include "fastq-reader.h"
#include "file.h"
#include "gzip-compressor.h"
#include "input-groups.h"
#include "lengths.h"
#include "logger.h"
#include "options.h"
#include "reads-reader.h"
#include "reference.h"
#include "scheme.h"
#include "strbuf.h"
#include "strutils.h"
#include "trimmer.h"

/* reads per thread in a batch */
#define BATCH_PER_THREAD 1000
#define COMPRESS_LEVEL 6
#define SUMMARY_FILENAME "summary.tsv"
#define LENGTHS_FILENAME "lengths.tsv"
/* the peaks of the lengths shown without references */
#define MAX_PEAKS 12

/* the paths of the outputs of the groups; an input which would be
   overwritten by an output is an error */
static char** output_paths(const input_group_t *groups, size_t n_groups,
    const char *outdir, Logger *logger) {
  char **paths;
  XMALLOC(paths, sizeof(char *) * (n_groups + 1));
  for (size_t i = 0; i < n_groups; i++) {
    paths[i] = input_group_output_path(outdir, groups[i].name);
    for (size_t j = 0; j < groups[i].n_files; j++)
      LCHECK(strcmp(groups[i].files[j], paths[i]) != 0, logger,
             "%s: the input would be overwritten by the output\n",
             groups[i].files[j]);
  }
  return paths;
}

static void delete_output_paths(char **paths, size_t n) {
  for (size_t i = 0; i < n; i++)
    free(paths[i]);
  free(paths);
}

/* the work of a thread on a part of a batch */
typedef struct {
  trimmer_t *trimmer;
  gzip_compressor_t *compressor;
  const fastq_batch_t *batch;
  size_t from, to;
  strbuf_t text, gz;
  counter_t *counts;
} worker_t;

static void* work(void *data) {
  worker_t *w = data;
  strbuf_clear(&w->text);
  strbuf_clear(&w->gz);
  trimmer_run(w->trimmer, w->batch, w->from, w->to, &w->text, w->counts);
  if (w->text.len > 0)
    gzip_compressor_run(w->compressor, w->text.s, w->text.len, &w->gz);
  return NULL;
}

/* trim the reads of <file>, appending them to <out> */
static void trim_file(const char *file, worker_t *workers, size_t n_workers,
    fastq_batch_t *batch, FILE *out, const char *out_path,
    Logger *logger) {
  reads_reader_t *reader = reads_reader_new(file, logger);
  pthread_t *threads;
  XMALLOC(threads, sizeof(pthread_t) * n_workers);
  size_t n;
  while ((n = reads_reader_read_batch(reader, batch,
                                      BATCH_PER_THREAD * n_workers)) > 0) {
    size_t per = (n + n_workers - 1) / n_workers;
    size_t used = 0;
    for (size_t t = 0; t < n_workers && t * per < n; t++, used++) {
      workers[t].batch = batch;
      workers[t].from = t * per;
      workers[t].to = (t + 1) * per < n ? (t + 1) * per : n;
    }
    if (used == 1)
      work(workers);
    else {
      for (size_t t = 0; t < used; t++)
        LCHECK(pthread_create(threads + t, NULL, work, workers + t) == 0,
               logger, "Cannot start a thread\n");
      for (size_t t = 0; t < used; t++)
        pthread_join(threads[t], NULL);
    }
    for (size_t t = 0; t < used; t++)
      LCHECK(fwrite(workers[t].gz.s, 1, workers[t].gz.len, out) ==
             workers[t].gz.len, logger, "%s: write error\n", out_path);
  }
  free(threads);
  reads_reader_delete(reader);
}

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
static void log_stats(const char *name, const counter_reads_t *s, bool refs,
    Logger *logger) {
  lprintf(logger, LOG_SUCCESS, "%s: %zu of %zu reads written (%zu with "
          "both primers, %zu with one, %zu with none; %zu too short, of "
          "which %zu dimers)\n", name, s->reads_out, s->reads_in,
          s->n_ends[2], s->n_ends[1], s->n_ends[0], s->too_short,
          s->dimers);
  if (refs)
    lprintf(logger, LOG_SUCCESS, "%s: ends started by internal priming: "
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
static size_t report_lengths(const lengths_t *l, const fasta_t *refs,
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
      lprintf(logger, LOG_SUCCESS, "%s: complete amplicons, lengths "
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
    lprintf(logger, LOG_SUCCESS, "%s: %s: %zu of %zu complete amplicons "
            "shorter than %.0f%% of the reference (%zu bases), most at "
            "%zu-%zu: internal priming or deletions\n", path, name,
            s.short_count, s.total, LENGTHS_SHORT_FRACTION * 100,
            refs->records[g].seq_len, s.short_bin * LENGTHS_BIN,
            s.short_bin * LENGTHS_BIN + LENGTHS_BIN - 1);
  }
  return n_short;
}

int main(int argc, char *argv[]) {
  options_t opts;
  Logger *logger;
  options_parse(&opts, argc, argv, &logger);
  if (opts.show_banner)
    banner_print(stderr);
  size_t n_groups;
  input_group_t *groups = input_groups_new(opts.inputs, opts.n_inputs,
                                           &n_groups, false, true, logger);
  char **paths = output_paths(groups, n_groups, opts.outdir, logger);
  scheme_t *scheme = scheme_new_builtin(opts.scheme, logger);
  fasta_t *refs = opts.references == NULL ? NULL :
                  fasta_read(opts.references, logger);
  reference_index_t *index = NULL;
  if (refs != NULL) {
    lprintf(logger, LOG_INFO, "%zu reference sequences read from %s\n",
            refs->n_records, opts.references);
    index = reference_index_new(refs);
  }
  trimmer_config_t config = {.scheme = scheme,
      .max_errors = opts.max_errors, .policy = opts.policy,
      .window = opts.window, .references = refs, .index = index,
      .exhaustive = opts.exhaustive_search,
      .probe_len = opts.probe_len,
      .probe_error_rate = opts.probe_error_rate,
      .terminal_window = opts.terminal_window,
      .min_length = opts.min_length, .tags = opts.tags};
  worker_t *workers;
  XMALLOC(workers, sizeof(worker_t) * opts.threads);
  for (size_t t = 0; t < opts.threads; t++) {
    workers[t].trimmer = trimmer_new(&config);
    workers[t].compressor = gzip_compressor_new(COMPRESS_LEVEL);
    strbuf_init(&workers[t].text);
    strbuf_init(&workers[t].gz);
    workers[t].counts = counter_new(scheme, refs);
  }
  /* the counts of each group, and of all reads */
  counter_t **group_counts;
  XMALLOC(group_counts, sizeof(counter_t *) * (n_groups + 1));
  counter_t *total = counter_new(scheme, refs);
  xmakedirs(opts.outdir, logger);
  char *summary_path = path_join(opts.outdir, SUMMARY_FILENAME);
  FILE *summary = xfopen(summary_path, NULL, "w", "summary table", logger);
  summary_header(summary);
  char *lengths_path = path_join(opts.outdir, LENGTHS_FILENAME);
  FILE *lengths_file = xfopen(lengths_path, NULL, "w", "lengths table",
                              logger);
  fputs("output\tgroup\tlength\treads\n", lengths_file);
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  for (size_t i = 0; i < n_groups; i++) {
    const input_group_t *g = groups + i;
    for (size_t t = 0; t < opts.threads; t++)
      counter_clear(workers[t].counts);
    FILE *out = xfopen(paths[i], NULL, "wb", "trimmed reads", logger);
    for (size_t j = 0; j < g->n_files; j++) {
      lprintf(logger, LOG_PHASE, "Trimming %s...\n", g->files[j]);
      trim_file(g->files[j], workers, opts.threads, &batch, out, paths[i],
                logger);
    }
    xfclose(&out, logger);
    counter_t *counts = group_counts[i] = counter_new(scheme, refs);
    for (size_t t = 0; t < opts.threads; t++)
      counter_merge(counts, workers[t].counts);
    counter_merge(total, counts);
    const counter_reads_t *stats = counter_reads(counts);
    log_stats(paths[i], stats, refs != NULL, logger);
    size_t n_short = report_lengths(counter_lengths(counts), refs, g->name,
                                    paths[i], lengths_file, logger);
    summary_line(summary, g->name, stats, refs != NULL, n_short);
  }
  xfclose(&summary, logger);
  xfclose(&lengths_file, logger);
  free(lengths_path);
  for (size_t i = 0; i < n_groups; i++)
    counter_delete(group_counts[i]);
  free(group_counts);
  counter_delete(total);
  fastq_batch_free(&batch);
  for (size_t t = 0; t < opts.threads; t++) {
    trimmer_delete(workers[t].trimmer);
    gzip_compressor_delete(workers[t].compressor);
    strbuf_free(&workers[t].text);
    strbuf_free(&workers[t].gz);
    counter_delete(workers[t].counts);
  }
  free(workers);
  free(summary_path);
  reference_index_delete(index);
  fasta_delete(refs);
  scheme_delete(scheme);
  delete_output_paths(paths, n_groups);
  input_groups_delete(groups, n_groups);
  delete_logger(logger);
  return EXIT_SUCCESS;
}
