/* (c) 2026 Giorgio Gonnella */

/*
   The short deletions in the flanks of a read with an internal deletion,
   next to its junction.

   Many reads of deletion-containing genomes carry, besides their
   junction, a second, short deletion (from a few to about 50 bases) a
   few dozen to about 200 bases from it. It is too short to be told from
   the indels of a read by the search of the junction (junction.h), so it
   stays in the flank. Where such reads are many, their flanks outvote
   the full templates in an assembly, and the deletion enters the
   consensus, where it can break the reading frame.

   The bases of each flank next to the junction (FLANK_DELETIONS_WINDOW,
   at most) are aligned to the reference, from the junction outwards, and
   each deletion of at least FLANK_DELETIONS_MIN bases found there is a
   cut point: the flank cut there lacks no base inside it, and the cut
   loses no base of the read. A difference with the reference which the
   sample itself has (e.g. a deletion in its neuraminidase stalk) only
   cuts the flanks there too: the full templates still carry it. Only
   the bases next to the junction are aligned, not the whole flank.
*/

#ifndef FLANK_DELETIONS_H
#define FLANK_DELETIONS_H

#include <stddef.h>
#include "junction.h"
#include "reference.h"

/* the bases of each flank aligned, from the junction */
#define FLANK_DELETIONS_WINDOW 250

/* the shortest deletion cut: the deletions of a read of Nanopore are
   rarely longer (a deletion which is only an error of the read is cut
   as well, which loses no base) */
#define FLANK_DELETIONS_MIN 6

/* the most cut points of a read */
#define FLANK_DELETIONS_MAX_CUTS 8

/*
   The cut points of the template <seq> (length <len>, in the orientation
   of the plus strand of the reference <ref_index> of <refs>) whose
   junction <j> was located (junction_locate()): the positions in the
   template of the first base after each short deletion of the flanks
   next to the junction, in increasing order, into <cuts> (at most
   FLANK_DELETIONS_MAX_CUTS); their number. None if the junction was not
   located, or where the bases differ too much from the reference.
*/
size_t flank_deletions_find(const reference_set_t *refs, const char *seq,
    size_t len, size_t ref_index, const junction_t *j, long *cuts);

/* Unit tests */
void flank_deletions_test(void);

#endif
