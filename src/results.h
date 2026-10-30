/* (c) 2026 Giorgio Gonnella */

/*
   The results of a run: the summary table (summary.tsv: the counts of
   each output), the lengths of the complete amplicons (lengths.tsv) with
   their check (see lengths.h), the messages of the counts, and the
   reports (run-report.h). Written at the end of a run, and in the live
   mode after each file (then without the messages).
*/

#ifndef RESULTS_H
#define RESULTS_H

#include <stdbool.h>
#include <stddef.h>
#include "counter.h"
#include "fasta.h"
#include "logger.h"
#include "run-report.h"

#define RESULTS_SUMMARY_FILE "summary.tsv"
#define RESULTS_LENGTHS_FILE "lengths.tsv"

/* a group of input files and its output */
typedef struct {
  const char *name;
  const char *output;        /* the path of the output file */
  const counter_t *counts;
  size_t n_files;
} results_group_t;

/* Log the counts <s> of the output <name> (<refs>: with references) */
void results_log_counts(const char *name, const counter_reads_t *s,
    bool refs, Logger *logger);

/*
   Write the tables and the reports of the <n> <groups> in <outdir>, the
   reports with the <total> counts; <running>: a run in progress (see
   run_report_write(); no messages)
*/
void results_write(const char *outdir, const results_group_t *groups,
    size_t n, const counter_t *total, const fasta_t *refs,
    const run_report_t *report, bool running, Logger *logger);

#endif
