/* (c) 2026 Giorgio Gonnella */

/*
   The plots of the HTML report, as inline SVG. Of a segment:
   - the coverage by read origin: a stacked profile along the reference,
     by category of the reads covering each position (read-cover.h; the
     strands together), where
     internal priming shows as a step up and internal deletions as a band
     over the two ends only;
   - the lengths of the complete amplicons, with the length of the
     reference marked: the amplicons shorter than it (internal priming,
     deletions) stand apart from those of the whole segment.
   Of a run: the median depth of each segment in each sample.
   The colors are those of the classes c-* of html-report-style.html.
*/

#ifndef HTML_REPORT_PLOTS_H
#define HTML_REPORT_PLOTS_H

#include <stdbool.h>
#include <stddef.h>
#include "counter.h"
#include "html-writer.h"
#include "lengths.h"

/* The legend of the coverage plots, drawn with <min_depth> */
void plots_coverage_legend(html_writer_t *w, size_t min_depth);

/* The coverage of the reference <ref>, of length <len>, by the located
   reads of <c>; <min_depth>: a line is drawn at this depth, and the
   stretches below it are marked under the axis (0: none); <mark_site>:
   the internal priming site with most read ends is marked */
void plots_coverage(html_writer_t *w, const counter_t *c, size_t ref,
    size_t len, size_t min_depth, bool mark_site);

/* The lengths of the complete amplicons of the group <group> of <l>;
   <ref_len>: the length of their reference, marked (0: none, the plot
   then ends at the longest amplicons) */
void plots_lengths(html_writer_t *w, const lengths_t *l, size_t group,
    size_t ref_len);

/* The median depths of the eight segments in <n> samples (<depths>: by
   sample, then by segment; <names>: of the samples), a dot for each
   sample over each segment, in the order of the samples, on a
   logarithmic axis, with a line at <min_depth>; the dot of a sample
   whose segments cannot be judged (<judged> false) is hollow */
void plots_segment_depths(html_writer_t *w, const size_t *depths,
    const char *const *names, const bool *judged, size_t n,
    size_t min_depth);

#endif
