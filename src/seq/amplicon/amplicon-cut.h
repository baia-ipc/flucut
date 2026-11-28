/* (c) 2026 Giorgio Gonnella */

/*
   The classification of the ends of a read against the references, and
   the part of the read to keep (see amplicon.h): a step of the analysis
   of a read (amplicon-analyzer.h).

   The probe of each end, read from its primer into the template, is
   located in the references: at a terminus of a segment the end is
   terminal, at an interior position it was internally primed. The part
   to keep follows from the primers found, the policy and the origin of
   each end: the conserved terminus is kept only at a terminal end of a
   segment the primer binds, with TRIM_TERMINUS; a primer dimer keeps
   nothing; a primer of the opposite orientation facing the layout is
   removed with what lies outside it.
*/

#ifndef AMPLICON_CUT_H
#define AMPLICON_CUT_H

#include <stddef.h>
#include "amplicon.h"

/* Classify the ends of the amplicon <a> found in <seq> (length <len>)
   against the references of the analyzer */
void amplicon_classify_ends(amplicon_analyzer_t *self, const char *seq,
    size_t len, amplicon_t *a);

/* The part of the read (length <len>) to keep, <from> to <to>, in the
   orientation in which the layout <a> was found, with the policy
   <policy> */
void amplicon_compute_keep(const amplicon_t *a, trim_policy_t policy,
    size_t len, long *from, long *to);

/* Limit <from> and <to> (see amplicon_compute_keep()) so that a primer
   of the layout <other>, of the opposite orientation, at an end of the
   read without a primer in <a>, is removed whole */
void amplicon_cut_facing_primer(const amplicon_t *a,
    const amplicon_t *other, size_t len, long *from, long *to);

#endif
