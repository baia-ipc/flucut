/* (c) 2026 Giorgio Gonnella */

/*
   The reads themselves in the HTML report: their number, bases, lengths
   (mean, median, N50, range) and qualities, before and after the
   trimming, as a table and as two plots of their distributions
   (read-stats.h); and what their titles say about the sequencing run
   (read-metadata.h).
*/

#ifndef HTML_REPORT_READS_H
#define HTML_REPORT_READS_H

#include "counter.h"
#include "html-writer.h"

/* The lengths and qualities of the reads of <c>, and the sequencing run */
void reads_write(html_writer_t *w, const counter_t *c);

#endif
