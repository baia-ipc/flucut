/* (c) 2026 Giorgio Gonnella */

/*
   The segments of a sample in the HTML report: one strip for each
   segment, with its main numbers (reference, depth, how much of it is
   well covered, the shares of its depth from internally primed reads and
   from reads with an internal deletion), its coverage by read origin and
   the lengths of its complete amplicons (html-report-plots.h). A
   reference whose description names no segment has its own strip, if it
   has reads.
*/

#ifndef HTML_REPORT_SEGMENTS_H
#define HTML_REPORT_SEGMENTS_H

#include "counter.h"
#include "html-writer.h"
#include "qc-summary.h"
#include "report.h"

/* The segments of the sample with the counts <c> and the summary <qc>
   (with references) */
void segments_write(html_writer_t *w, const report_t *r, const counter_t *c,
    const qc_summary_t *qc);

/* The segment <segment> (1 to 8) in all the samples of the run: the
   strip of the segment for each sample, named after the sample and
   linked to its part of the report (#s<i>) */
void segments_write_across(html_writer_t *w, const report_t *r,
    int segment);

/* How the bases of the located reads divide among the segments of the
   summary <qc>: a column for each segment (nothing if no reads were
   located) */
void segments_write_balance(html_writer_t *w, const qc_summary_t *qc);

/* Without references: the lengths of all the complete amplicons */
void segments_write_lengths(html_writer_t *w, const counter_t *c);

#endif
