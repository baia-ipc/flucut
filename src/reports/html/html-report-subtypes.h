/* (c) 2026 Giorgio Gonnella */

/*
   The reads of each HA and NA subtype of a sample in the HTML report
   (subtype-reads.h): for the HA and the NA, the subtype of most reads
   and the second one, with their reads, and, when the reads of the
   second one are above the threshold, a remark which says what was
   observed and what can cause it.
*/

#ifndef HTML_REPORT_SUBTYPES_H
#define HTML_REPORT_SUBTYPES_H

#include "counter.h"
#include "html-writer.h"

/* The card of the sample with the counts <c> (nothing if there is
   nothing to say, e.g. influenza B or the counts of several samples) */
void subtypes_write(html_writer_t *w, const counter_t *c);

#endif
