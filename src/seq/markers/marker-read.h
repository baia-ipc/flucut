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

   <cuts>: the <n_cuts> positions at which the read is not contiguous
   on its reference - the junction of an internal deletion and the short
   deletions beside it (amplicon.h), each the first base after the
   missing stretch. A marker is read from such a read only where its
   anchor and its window lie between two of them. Without this, the
   sites which happen to lie beside the junction of a defective genome
   are read across it, and a sample full of such genomes reports the
   residue of another part of the segment.

   The cuts are those of the read as it arrived, whether or not the
   output is written in pieces (--deletions): what a marker says is
   about the sample, not about how the reads are written.
*/
void marker_read_sites(const marker_reference_t *const *chosen, long ref,
    const long *cuts, size_t n_cuts, const char *read, size_t len,
    double max_error_rate, marker_counts_t *counts);

#endif
