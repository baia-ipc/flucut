/* (c) 2026 Giorgio Gonnella */

/*
   A sample in the HTML report (a group of input files, or all the reads
   of a run of one sample): its state and what deserves attention in it,
   its main numbers as tiles, how its reads divide (by primers found, by
   origin), its internal priming sites, and its segments
   (html-report-segments.h).
*/

#ifndef HTML_REPORT_SAMPLE_H
#define HTML_REPORT_SAMPLE_H

#include "counter.h"
#include "html-report-markup.h"
#include "html-writer.h"
#include "qc-summary.h"
#include "report.h"

/* The badge of the state of a sample: its segments complete, of those
   expected (<brief>: as "3/8", otherwise "3/8 segments complete"). It
   is neutral for a sample not complete whose segments cannot be judged
   (few reads, reads not located), or not yet (<r>: a run in progress
   which is still adding reads to its samples, in the live mode) */
void sample_badge(html_writer_t *w, const report_t *r,
    const qc_summary_t *qc, bool brief);

/* The subtype of the sample summarized as <qc>, followed by a question
   mark if its HA and NA do not support it well (one of them not found or
   at a low depth, or a sample whose segments cannot be judged) */
typedef struct {
  char s[QC_SUBTYPE_SIZE + 2];
} sample_subtype_t;
sample_subtype_t sample_subtype(const qc_summary_t *qc);

/* A number of bases in a short form, e.g. "11.6 Mb" */
typedef struct {
  char s[24];
} sample_bases_t;
sample_bases_t sample_bases(size_t bases);

/* The tile of the lengths of the reads of the input of <c>: the median,
   with the N50 and the median quality */
void sample_length_tile(html_writer_t *w, const counter_t *c);

/* The line under the name of a sample with the counts <c>: its state,
   its subtype and what deserves attention (nothing without references
   and without anything to say) */
void sample_status_line(html_writer_t *w, const report_t *r,
    const counter_t *c);

/* The content of a sample with the counts <c> */
void sample_write(html_writer_t *w, const report_t *r, const counter_t *c);

#endif
