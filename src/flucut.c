/* (c) 2026 Giorgio Gonnella */

/*
   flucut: removal of the primers of the influenza A multisegment RT-PCR
   from the reads, before the assembly (see README.md).

   Each input (a file of reads, or a directory of files of reads, e.g. a
   barcode directory) gives one gzipped FASTQ file in the output
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
#include "fasta.h"
#include "fastq-reader.h"
#include "file.h"
#include "gzip-compressor.h"
#include "logger.h"
#include "options.h"
#include "reads-reader.h"
#include "scheme.h"
#include "strbuf.h"
#include "strutils.h"
#include "trimmer.h"

/* reads per thread in a batch */
#define BATCH_PER_THREAD 1000
#define COMPRESS_LEVEL 6
#define SUMMARY_FILENAME "summary.tsv"

/* an output: its name, file and inputs */
typedef struct {
  char *name;
  char *path;
  char **files;
  size_t n_files;
} output_t;

/* list_dir() filter: the files of reads */
static bool is_reads_file(const char *name, const char *path, void *data) {
  (void)data;
  return reads_reader_is_reads_filename(name) && !is_directory(path);
}

/* the outputs of the inputs; two inputs with the same output name are an
   error, as one would overwrite the other */
static output_t* plan_outputs(const options_t *opts, Logger *logger) {
  output_t *outputs;
  XMALLOC(outputs, sizeof(output_t) * opts->n_inputs);
  for (size_t i = 0; i < opts->n_inputs; i++) {
    const char *in = opts->inputs[i];
    output_t *o = outputs + i;
    o->name = path_basename(in);
    if (is_directory(in)) {
      o->files = list_dir(in, is_reads_file, NULL, &o->n_files, logger);
      LCHECK(o->n_files > 0, logger, "%s: no file of reads in the "
             "directory (.fastq, .fq, optionally .gz, or .bam)\n", in);
    } else {
      LCHECK(reads_reader_is_reads_filename(o->name), logger,
             "%s: not a file of reads (.fastq, .fq, optionally .gz, or "
             ".bam)\n", in);
      o->name[strlen(o->name) - reads_reader_filename_suffix_len(o->name)]
        = '\0';
      o->n_files = 1;
      XMALLOC(o->files, sizeof(char *));
      o->files[0] = xstrdup(in);
    }
    for (size_t j = 0; j < i; j++)
      LCHECK(strcmp(outputs[j].name, o->name) != 0, logger,
             "%s and %s: both would be written to %s.fastq.gz\n",
             opts->inputs[j], in, o->name);
    size_t len = strlen(o->name) + sizeof(".fastq.gz");
    char *filename;
    XMALLOC(filename, len);
    snprintf(filename, len, "%s.fastq.gz", o->name);
    o->path = path_join(opts->outdir, filename);
    free(filename);
    for (size_t j = 0; j < o->n_files; j++)
      LCHECK(strcmp(o->files[j], o->path) != 0, logger,
             "%s: the input would be overwritten by the output\n",
             o->files[j]);
  }
  return outputs;
}

static void delete_outputs(output_t *outputs, size_t n) {
  for (size_t i = 0; i < n; i++) {
    free(outputs[i].name);
    free(outputs[i].path);
    for (size_t j = 0; j < outputs[i].n_files; j++)
      free(outputs[i].files[j]);
    free(outputs[i].files);
  }
  free(outputs);
}

/* the work of a thread on a part of a batch */
typedef struct {
  trimmer_t *trimmer;
  gzip_compressor_t *compressor;
  const fastq_batch_t *batch;
  size_t from, to;
  strbuf_t text, gz;
  trimmer_stats_t stats;
} worker_t;

static void* work(void *data) {
  worker_t *w = data;
  strbuf_clear(&w->text);
  strbuf_clear(&w->gz);
  trimmer_run(w->trimmer, w->batch, w->from, w->to, &w->text, &w->stats);
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
        "left_unknown\tright_terminal\tright_internal\tright_unknown\n",
        f);
}

static void summary_line(FILE *f, const char *name,
    const trimmer_stats_t *s, bool refs) {
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
  fputs("\n", f);
}

/* log the counts of an output */
static void log_stats(const char *name, const trimmer_stats_t *s, bool refs,
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

int main(int argc, char *argv[]) {
  options_t opts;
  Logger *logger;
  options_parse(&opts, argc, argv, &logger);
  output_t *outputs = plan_outputs(&opts, logger);
  scheme_t *scheme = scheme_new_builtin(opts.scheme, logger);
  fasta_t *refs = opts.references == NULL ? NULL :
                  fasta_read(opts.references, logger);
  if (refs != NULL)
    lprintf(logger, LOG_INFO, "%zu reference sequences read from %s\n",
            refs->n_records, opts.references);
  trimmer_config_t config = {.scheme = scheme,
      .max_errors = opts.max_errors, .policy = opts.policy,
      .window = opts.window, .references = refs,
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
  }
  xmakedirs(opts.outdir, logger);
  char *summary_path = path_join(opts.outdir, SUMMARY_FILENAME);
  FILE *summary = xfopen(summary_path, NULL, "w", "summary table", logger);
  summary_header(summary);
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  for (size_t i = 0; i < opts.n_inputs; i++) {
    output_t *o = outputs + i;
    for (size_t t = 0; t < opts.threads; t++)
      memset(&workers[t].stats, 0, sizeof(trimmer_stats_t));
    FILE *out = xfopen(o->path, NULL, "wb", "trimmed reads", logger);
    for (size_t j = 0; j < o->n_files; j++) {
      lprintf(logger, LOG_PHASE, "Trimming %s...\n", o->files[j]);
      trim_file(o->files[j], workers, opts.threads, &batch, out, o->path,
                logger);
    }
    xfclose(&out, logger);
    trimmer_stats_t stats;
    memset(&stats, 0, sizeof(stats));
    for (size_t t = 0; t < opts.threads; t++)
      trimmer_stats_add(&stats, &workers[t].stats);
    summary_line(summary, o->name, &stats, refs != NULL);
    log_stats(o->path, &stats, refs != NULL, logger);
  }
  xfclose(&summary, logger);
  fastq_batch_free(&batch);
  for (size_t t = 0; t < opts.threads; t++) {
    trimmer_delete(workers[t].trimmer);
    gzip_compressor_delete(workers[t].compressor);
    strbuf_free(&workers[t].text);
    strbuf_free(&workers[t].gz);
  }
  free(workers);
  free(summary_path);
  fasta_delete(refs);
  scheme_delete(scheme);
  delete_outputs(outputs, opts.n_inputs);
  delete_logger(logger);
  return EXIT_SUCCESS;
}
