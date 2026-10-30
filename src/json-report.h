/* (c) 2026 Giorgio Gonnella */

/*
   JSON summary (from Lumipore): the statistics of a run (the same as in
   the HTML report, see report.h), in a machine-readable form, e.g. for
   pipelines. The format is described in docs/json_summary.md. The file
   is written to a temporary file and then renamed, so that it is never
   read partially written.
*/

#ifndef JSON_REPORT_H
#define JSON_REPORT_H

#include "logger.h"
#include "report.h"

/* version of the format of the JSON summary */
#define JSON_SUMMARY_FORMAT 1

/* maximum number of internal priming sites output for each group */
#define JSON_MAX_SITES 1000

/* Write the summary to <path>; exits with an error message on failure */
void json_report_write(const char *path, const report_t *report,
    Logger *logger);

#endif
