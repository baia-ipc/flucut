/* (c) 2026 Giorgio Gonnella */

/*
   The location on the references of what the probes of the ends of a
   read did not locate (see amplicon.h: "Placement", "Ends located from
   inside", "Internal deletions"): a step of the analysis of a read
   (amplicon-analyzer.h), after the part to keep is known. None of it
   changes the trimming: it is reported, for the coverage.

   - a read of which no end was located is placed by a probe from inside
     its kept part (placement.h);
   - an end not located by its own probe is located by probes from
     further inside the read, up to AMPLICON_INNER_PROBES of them;
   - a read whose two ends map to the same reference, farther apart than
     the template between them is long, has the junction of its internal
     deletion located (junction.h).

   Each function has no effect without references.
*/

#ifndef AMPLICON_LOCATE_H
#define AMPLICON_LOCATE_H

#include <stddef.h>
#include "amplicon.h"

/* Place on the references the kept part of the read <seq> analyzed in
   <a>, if no end of it was located */
void amplicon_place(amplicon_analyzer_t *self, const char *seq,
    amplicon_t *a);

/* Locate from inside the read <seq> (length <len>) the ends of <a> which
   were not located */
void amplicon_locate_inner_ends(amplicon_analyzer_t *self, const char *seq,
    size_t len, amplicon_t *a);

/* Locate the junction of the read <seq> (length <len>) analyzed in <a>,
   if its two ends map to the same reference farther apart than the
   template between them is long */
void amplicon_find_junction(amplicon_analyzer_t *self, const char *seq,
    size_t len, amplicon_t *a);

#endif
