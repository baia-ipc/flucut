/* (c) 2026 Giorgio Gonnella */

/*
   The run in the head of the HTML report: its main settings on one line,
   its state while it is in progress, and all the settings (version,
   command line, scheme, references, options) in a table shown when the
   reader opens it.
*/

#ifndef HTML_REPORT_RUN_H
#define HTML_REPORT_RUN_H

#include "html-writer.h"
#include "report.h"

void run_write(html_writer_t *w, const report_t *r);

#endif
