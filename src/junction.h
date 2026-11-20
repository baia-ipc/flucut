/* (c) 2026 Giorgio Gonnella */

/*
   The junction of a read with an internal deletion.

   The two ends of such a read (e.g. of a defective interfering genome)
   map to the two ends of a segment, much farther apart than the read is
   long: the read is two flanks of the segment joined at a junction, and
   what lies between them on the reference is missing. Where the junction
   is cannot be told from the ends; it is located here, so that each
   flank is counted where it lies (see counter.h).

   Probes from inside the read are located in the reference of the read
   (reference_set_locate_in()): a probe of the left flank maps where the
   left end of the read puts it, a probe of the right flank where the
   right end puts it, and the junction is between the last probe of the
   left flank and the first of the right one. It is searched by
   bisection, so that the number of probes grows with the logarithm of
   the length of the read. Between the two probes, the bases are aligned
   to the reference from both sides, along the two flanks: the junction
   is where the bases before it and the bases after it differ the least
   from the reference. If a probe does not tell its flank (it contains
   the junction, or the sample differs too much from the reference
   there), the bisection ends there, and the bases between the probes
   known so far are aligned, at most the whole read.

   Only the probes and the bases next to the junction are aligned, not
   the read. Bases which are in the reference on both sides of the
   deletion (a short repeat at the junction) belong to either flank: they
   are given to the right one. A difference between the sample and the
   reference next to the junction (bases which the sample lacks) cannot
   be told from the deletion, and is part of it.
*/

#ifndef JUNCTION_H
#define JUNCTION_H

#include <stdbool.h>
#include <stddef.h>
#include "reference.h"

/* a deletion smaller than this is not told from the indels of a read */
#define JUNCTION_MIN_DELETION 100

/* the junction of a read on its reference */
typedef struct {
  bool found;
  /* the flanks on the plus strand of the reference: the left one ends
     before <left_to>, the right one starts at <right_from>; the bases
     [left_to, right_from) of the reference are missing in the read */
  long left_to, right_from;
} junction_t;

/*
   Whether a template of <len> bases whose ends mapped <distance> bases
   apart has an internal deletion: it is shorter than the distance by
   more than JUNCTION_MIN_DELETION and a tenth of the distance (the
   allowance for the indels of a long read)
*/
bool junction_expected(long len, long distance);

/*
   Locate the junction of the template <seq> (length <len>), given in the
   orientation of the plus strand of the reference <ref_index> of <refs>,
   whose first base is at <first> and whose last base is at <last> of the
   reference; into <out> (found: false if the template has no deletion,
   see junction_expected(), or its bases do not fit the two flanks)
*/
void junction_locate(reference_set_t *refs, const char *seq, size_t len,
    size_t ref_index, long first, long last, junction_t *out);

/*
   The differences (edit distance) between the sequence <seq> (length
   <len>) and the reference <ref_index> of <refs>, aligned from one of
   its ends: forwards from its first base at <anchor> of the reference,
   or (<backwards>) from its last base at <anchor> towards the start of
   both; <seq> is given in the orientation of the plus strand. For
   telling whether bases of a read continue the reference at a position
*/
long junction_flank_differences(reference_set_t *refs, const char *seq,
    size_t len, size_t ref_index, long anchor, bool backwards);

/* Unit tests */
void junction_test(void);

#endif
