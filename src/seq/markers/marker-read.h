/* (c) 2026 Giorgio Gonnella */

/*
   Reading the markers (markers.h) out of one read.

   Which reads are read is left to the caller, because it is not the
   same on the two platforms: a read of a whole amplicon (Nanopore) is
   located on a reference, so only the markers of that reference are
   looked for in it, while a read of a fragment (Illumina) which covers
   a marker usually carries no primer and is not located at all, so
   every marker is looked for and the anchors themselves say which ones
   the read covers.
*/

#ifndef MARKER_READ_H
#define MARKER_READ_H

#include <stddef.h>
#include "marker-counts.h"
#include "markers.h"

/*
   Read the markers of <chosen> (MARKER_N_MARKERS anchors, NULL where
   the sample has none, see marker_anchors_chosen()) out of the <len>
   bases of <read>, with at most <max_error_rate> of its length as
   errors in an anchor, and count them in <counts>.

   <ref>: the reference the read was located on, so that only the
   markers taken from it are looked for; negative where the read was
   not located, and then every marker is.

   <breaks>: where the read does not follow its reference. A marker is
   read only from a stretch which has none of them, because the bases
   there are not the ones the position is named after:

   - a cut is a position at which the read is not contiguous - the
     junction of an internal deletion and the short deletions beside it
     (amplicon.h). A read with one is two flanks of the segment joined,
     and the bases before the junction come from somewhere else;
   - a blocked stretch is one whose bases are not the template: a
     primer which matched inside the read (amplicon.h). A read joined
     at a segment terminus carries the primer there, and a window read
     across it reads the primer: the terminus AGCAAAAGCAGG gives K at a
     marker which reads one codon, and SKS at one which reads three.

   The breaks are those of the read as it arrived, whether or not the
   output is written in pieces (--deletions): what a marker says is
   about the sample, not about how the reads are written.
*/
/* a stretch, 0-based with the last position included */
typedef struct {
  long from, to;
} marker_span_t;

/* how far the stretch of the reference a read covers may be wrong
   before a marker is refused for lying outside it: the span is
   estimated from where the ends of the read mapped and from its length,
   not from an alignment (read-cover.h) */
#define MARKER_PLACEMENT_SLACK 100

typedef struct {
  const long *cuts;                /* in the read */
  size_t n_cuts;
  const marker_span_t *blocked;    /* in the read */
  size_t n_blocked;
  /* the stretches of the reference the read covers (one, or two for a
     read with an internal deletion), in the reference. The anchor must
     lie in one of them: a read whose anchor matched where the read does
     not reach is joined there, and the bases beside the match are of
     another part of the segment. None: the read was not located, and
     the anchor itself says where it is */
  const marker_span_t *covered;
  size_t n_covered;
} marker_breaks_t;

void marker_read_sites(const marker_reference_t *const *chosen, long ref,
    const marker_breaks_t *breaks, const char *read, size_t len,
    double max_error_rate, marker_counts_t *counts);

#endif
