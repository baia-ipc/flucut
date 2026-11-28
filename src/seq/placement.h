/* (c) 2026 Giorgio Gonnella */

/*
   Approximate placement of a read on the references, for the reads of
   which no end was located (no primer at either end, or a primer whose
   probe did not map, see amplicon.h): such a read is still a piece of a
   segment, and belongs to the coverage of the segment (see counter.h).

   A probe is taken from inside the read and located by its k-mers
   (reference_set_locate_seeded()); the span of the read on the reference
   follows from the position of the probe in the read and the length of
   the read. Up to PLACEMENT_N_PROBES probes are tried, from the middle of
   the read outwards, until one is located.

   This is not a mapping of the read: only the probe is aligned, and the
   span ignores the indels of the rest of the read (the probe in the
   middle halves their effect at the ends). The cost for a read is that
   of locating a few probes by their k-mers, and does not depend on the
   length of the read or on the number of references.
*/

#ifndef PLACEMENT_H
#define PLACEMENT_H

#include <stdbool.h>
#include <stddef.h>
#include "reference.h"

/* the greatest number of probes tried for a read */
#define PLACEMENT_N_PROBES 3

/* a read placed on a reference */
typedef struct {
  bool found;
  size_t ref_index;    /* the reference */
  long ref_len;        /* its length */
  /* the read is on the minus strand of the reference */
  bool reverse;
  /* the span of the read on the plus strand of the reference:
     [from, to), not limited to the reference (a read which is longer
     than the reference on one side has from < 0 or to > ref_len) */
  long from, to;
} placement_t;

/*
   Place the sequence <seq> (length <len>; the part of a read without its
   primers and adapters) on the references of <refs>, into <out>
   (found: false if no probe is located, or the sequence is shorter than
   REFERENCE_MIN_PROBE_LEN)
*/
void placement_locate(reference_set_t *refs, const char *seq, size_t len,
    placement_t *out);

/* Unit tests */
void placement_test(void);

#endif
