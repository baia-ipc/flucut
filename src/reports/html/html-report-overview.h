/* (c) 2026 Giorgio Gonnella */

/*
   The overview of a run of several samples in the HTML report: how many
   samples are in each state, the throughput of the run as tiles, and a
   table with one row for each sample: its state, subtype and reads, a
   cell for each segment (how much of it is well covered), the segment
   with most internal priming and that with most internal deletions. The
   names of the samples link to their parts of the report (#s<i>).
   Nothing but the throughput is pooled across the samples.
*/

#ifndef HTML_REPORT_OVERVIEW_H
#define HTML_REPORT_OVERVIEW_H

#include "html-writer.h"
#include "report.h"

void overview_write(html_writer_t *w, const report_t *r);

#endif
