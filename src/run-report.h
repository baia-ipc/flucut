/* (c) 2026 Giorgio Gonnella */

/*
   The reports of a run (from Lumipore): the context of the HTML report
   and of the JSON summary (command line, options, scheme, references),
   and their writing, at the end of the run (or, in a run in progress,
   after each input file).
*/

#ifndef RUN_REPORT_H
#define RUN_REPORT_H

#include <stdbool.h>
#include <stddef.h>
#include "counter.h"
#include "fasta.h"
#include "logger.h"
#include "options.h"
#include "report.h"
#include "scheme.h"

/* files in the output directory, reload interval of the HTML report in
   a run in progress (seconds) */
#define RUN_REPORT_HTML_FILE "flucut_report.html"
#define RUN_REPORT_JSON_FILE "flucut_summary.json"
#define RUN_REPORT_REFRESH 30

typedef struct run_report_t run_report_t;

/*
   Create the run report for the options <opts>, the <scheme> and the
   <references> (NULL: none), which must exist as long as the report is
   used; the command line (<argc>, <argv>) and the <version> are output
   in the reports
*/
run_report_t* run_report_new(const options_t *opts, const char *version,
    const scheme_t *scheme, const fasta_t *references, int argc,
    char *argv[], Logger *logger);

/* Set the errors in a primer decided automatically (<setting> must exist
   as long as the report is used; it can change, e.g. in a run in
   progress) */
void run_report_set_max_errors(run_report_t *self,
    const report_setting_t *setting);

/*
   Write the reports of the <n> <sections> (groups of input files) and
   of the <total> counts; <running>: in a run in progress, before its end
*/
void run_report_write(const run_report_t *self,
    const report_section_t *sections, size_t n, const counter_t *total,
    bool running, Logger *logger);

/* Delete the report (NULL is accepted) */
void run_report_delete(run_report_t *self);

#endif
