/* (c) 2026 Giorgio Gonnella */

/*
   flucut: removal of the primers of the influenza A multisegment RT-PCR
   from the reads, before the assembly (see README.md).

   Each group of input files (input-groups.h: a file of reads, a
   directory of files of reads, or each barcode directory of a MinKNOW
   fastq_pass directory) gives one gzipped FASTQ file in the output
   directory, and a line of the summary table (summary.tsv); the
   reports (HTML and JSON) are written at the end. The files are trimmed
   by worker threads (worker.h), whose output does not depend on their
   number.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "auto-settings.h"
#include "banner.h"
#include "counter.h"
#include "fasta.h"
#include "file.h"
#include "input-groups.h"
#include "logger.h"
#include "options.h"
#include "reference.h"
#include "results.h"
#include "run-report.h"
#include "scheme.h"
#include "trimmer.h"
#include "worker.h"


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

/* log the scheme and the input groups (verbose) */
static void log_setup(const scheme_t *scheme, const input_group_t *groups,
    size_t n_groups, Logger *logger) {
  lprintf(logger, LOG_INFO, "Scheme %s: left", scheme->name);
  for (size_t i = 0; i < scheme->n_left; i++)
    lnprintf(logger, LOG_INFO, " %s (%s)", scheme->left[i].name,
             scheme->left[i].seq);
  lnprintf(logger, LOG_INFO, "; right");
  for (size_t i = 0; i < scheme->n_right; i++)
    lnprintf(logger, LOG_INFO, " %s (%zu bases)", scheme->right[i].name,
             scheme->right[i].seq_len);
  lnprintf(logger, LOG_INFO, "\n");
  lprintf(logger, LOG_INFO, "Input groups: %zu\n", n_groups);
  for (size_t i = 0; i < n_groups; i++)
    lprintf(logger, LOG_INFO, "  %s: %zu file%s\n", groups[i].name,
            groups[i].n_files, groups[i].n_files == 1 ? "" : "s");
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
  log_setup(scheme, groups, n_groups, logger);
  fasta_t *refs = opts.references == NULL ? NULL :
                  fasta_read(opts.references, logger);
  reference_index_t *index = NULL;
  if (refs != NULL) {
    lprintf(logger, LOG_INFO, "%zu reference sequences read from %s\n",
            refs->n_records, opts.references);
    index = reference_index_new(refs);
  }
  auto_settings_t autos;
  auto_settings_init(&autos, &opts);
  trimmer_config_t config = {.scheme = scheme,
      .max_errors = autos.max_errors.value, .policy = opts.policy,
      .window = opts.window, .references = refs, .index = index,
      .exhaustive = opts.exhaustive_search,
      .probe_len = opts.probe_len,
      .probe_error_rate = opts.probe_error_rate,
      .terminal_window = opts.terminal_window,
      .adapters = opts.adapters, .min_length = opts.min_length,
      .tags = opts.tags};
  /* -e auto: from the first file */
  for (size_t i = 0; i < n_groups && auto_settings_pending(&autos); i++)
    if (groups[i].n_files > 0)
      auto_settings_decide(&autos, &config, groups[i].files[0], logger);
  worker_t *workers;
  XMALLOC(workers, sizeof(worker_t) * opts.threads);
  for (size_t t = 0; t < opts.threads; t++)
    worker_init(workers + t, &config);
  /* the counts of each group, and of all reads */
  counter_t **group_counts;
  XMALLOC(group_counts, sizeof(counter_t *) * (n_groups + 1));
  counter_t *total = counter_new(scheme, refs);
  run_report_t *report = run_report_new(&opts, VERSION, scheme, refs, argc,
                                        argv, logger);
  run_report_set_max_errors(report, &autos.max_errors);
  results_group_t *results;
  XMALLOC(results, sizeof(results_group_t) * (n_groups + 1));
  for (size_t i = 0; i < n_groups; i++) {
    const input_group_t *g = groups + i;
    counter_t *counts = group_counts[i] = counter_new(scheme, refs);
    FILE *out = xfopen(paths[i], NULL, "wb", "trimmed reads", logger);
    for (size_t j = 0; j < g->n_files; j++) {
      lprintf(logger, LOG_PHASE, "Trimming %s...\n", g->files[j]);
      worker_process_file(workers, opts.threads, g->files[j], out, paths[i],
                          counts, logger);
    }
    xfclose(&out, logger);
    counter_merge(total, counts);
    results[i] = (results_group_t){g->name, paths[i], counts, g->n_files};
  }
  results_write(opts.outdir, results, n_groups, total, refs, report, false,
                logger);
  free(results);
  run_report_delete(report);
  auto_settings_free(&autos);
  for (size_t i = 0; i < n_groups; i++)
    counter_delete(group_counts[i]);
  free(group_counts);
  counter_delete(total);
  for (size_t t = 0; t < opts.threads; t++)
    worker_free(workers + t);
  free(workers);
  reference_index_delete(index);
  fasta_delete(refs);
  scheme_delete(scheme);
  delete_output_paths(paths, n_groups);
  input_groups_delete(groups, n_groups);
  delete_logger(logger);
  return EXIT_SUCCESS;
}
