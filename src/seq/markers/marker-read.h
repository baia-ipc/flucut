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
typedef struct {
  long from, to;   /* 0-based, the last one included */
} marker_blocked_t;

typedef struct {
  const long *cuts;
  size_t n_cuts;
  const marker_blocked_t *blocked;
  size_t n_blocked;
} marker_breaks_t;

void marker_read_sites(const marker_reference_t *const *chosen, long ref,
    const marker_breaks_t *breaks, const char *read, size_t len,
    double max_error_rate, marker_counts_t *counts);

#endif
