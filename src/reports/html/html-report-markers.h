/* (c) 2026 Giorgio Gonnella */

/*
   The markers of a sample in the HTML report: one row per marker the
   reads gave (markers.h), with what the reads carry at it, and, where
   more than one value is carried by enough reads, the values with
   their shares, so that a residue present in a part of the reads only
   is visible beside the one of the majority.
*/

#ifndef HTML_REPORT_MARKERS_H
#define HTML_REPORT_MARKERS_H

#include "counter.h"
#include "html-page.h"

/* Write the card of the markers of <c>; nothing where none was read */
void markers_write(html_writer_t *w, const counter_t *c);

/* Whether the table has a marker on <segment>, so that the table of
   that segment across the samples has a column for them */
bool markers_on_segment(int segment);

/*
   Write the cell of the markers of <segment> in <c>: each one with the
   value the reads carry, and a mark where it deserves a look
*/
void markers_cell(html_writer_t *w, const counter_t *c, int segment);

#endif
