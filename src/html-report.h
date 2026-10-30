/* (c) 2026 Giorgio Gonnella */

/*
   HTML report (from Lumipore): a self-contained HTML file (no external
   resources, plots as inline SVG, light and dark themes), with the
   statistics of a run:
   - the run parameters;
   - for all reads and for each group of input files (e.g. a barcode):
     the reads at each step (primers found, dimers, too short, written),
     the primers found at each end, the internal priming sites (with
     references), the complete amplicons and their lengths.
   The content is described in report.h. The file is written to a
   temporary file and then renamed, so that a browser never shows a
   partial report.
*/

#ifndef HTML_REPORT_H
#define HTML_REPORT_H

#include "logger.h"
#include "report.h"

/* the internal priming sites shown for each group */
#define HTML_REPORT_MAX_SITES 20

/* Write the report to <path>; exits with an error message on failure */
void html_report_write(const char *path, const report_t *report,
    Logger *logger);

#endif
