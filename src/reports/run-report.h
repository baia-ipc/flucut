/* (c) 2026 Giorgio Gonnella */

/*
   The reports of a run: the context of the HTML pages
   and of the JSON summary (command line, options, scheme, references),
   and their writing. The dashboard and the JSON summary are written
   while the run is in progress (after each group of input files; in the
   live mode, after each file) and at its end; the HTML report at its
   end.
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

/* files in the output directory, reload interval of the dashboard in a
   run in progress (seconds) */
#define RUN_REPORT_HTML_FILE "flucut_report.html"
#define RUN_REPORT_DASHBOARD_FILE "flucut_dashboard.html"
#define RUN_REPORT_JSON_FILE "flucut_summary.json"
#define RUN_REPORT_REFRESH 30

typedef struct run_report_t run_report_t;

/*
   Create the run report for the options <opts>, the <scheme> and the
   <references> (NULL: none), which must exist as long as the report is
   used; the command line (<argc>, <argv>) and the <version> are output
   in the reports
*/
/* The options of the run */
const options_t* run_report_options(const run_report_t *self);

run_report_t* run_report_new(const options_t *opts, const char *version,
    const scheme_t *scheme, const fasta_t *references, int argc,
    char *argv[], Logger *logger);

/* Set the errors in a primer decided automatically (<setting> must exist
   as long as the report is used; it can change, e.g. in a run in
   progress) */
void run_report_set_max_errors(run_report_t *self,
    const report_setting_t *setting);

/* Set the reload interval of the dashboard while a run is in progress
   (live mode; seconds); see live.c for the adaptive rule */
void run_report_set_refresh(run_report_t *self, long refresh);

/*
   Write the reports of the <n> <sections> (groups of input files) and
   of the <total> counts; <running>: in a run in progress, before its end.
   The resources used by the run (run-resources.h) are measured at each
   writing.
*/
void run_report_write(run_report_t *self,
    const report_section_t *sections, size_t n, const counter_t *total,
    bool running, Logger *logger);

/* Tell the user where the dashboard is and how to follow the run with it
   (to be called once it has been written) */
void run_report_announce_dashboard(const run_report_t *self,
    Logger *logger);

/* Delete the report (NULL is accepted) */
void run_report_delete(run_report_t *self);

#endif
