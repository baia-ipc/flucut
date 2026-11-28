/* (c) 2026 Giorgio Gonnella */

/*
   The search of the primers of the scheme at the ends of a read, and of
   the primer matches inside it (see amplicon.h): a step of the analysis
   of a read (amplicon-analyzer.h).

   The primer of an end is searched only within the window at that end of
   the read; of the primers of the end, the one which reaches furthest
   out is taken, then the one with fewer errors, then the longer one. At
   an end, the outer bases of a primer may be missing, after a flank or
   beyond the end of the read (amplicon_missing_costs()); the primer may
   overhang the read only if it is not found inside it. The matches
   between the primers of the two ends are collected separately
   (internal-matches.h).
*/

#ifndef AMPLICON_ENDS_H
#define AMPLICON_ENDS_H

#include <stdbool.h>
#include <stddef.h>
#include "amplicon.h"

/* The errors allowed in a primer of length <seq_len>: <max_errors>
   bases, or that fraction of its length if smaller than 1 */
long amplicon_max_errors_of(double max_errors, size_t seq_len);

/* The costs <skip> and <edge> (<len> + 1 values each) of the missing
   outer bases of a primer of length <len>, with <k> errors allowed in
   the whole primer (see amplicon-ends.c) */
void amplicon_missing_costs(double max_errors, size_t len, long k,
    long *skip, long *edge);

/* Search the primers of the left end at the start of <seq> (length
   <len>) and those of the right end at its end, within the windows
   <left_window> and <right_window>, into the ends of <out> */
void amplicon_search_layout(amplicon_analyzer_t *self, const char *seq,
    size_t len, long left_window, long right_window, amplicon_t *out);

/* Collect into <out> the matches of the primers of the scheme which lie
   between the primers of the two ends of <out> */
void amplicon_search_internal(amplicon_analyzer_t *self, const char *seq,
    size_t len, amplicon_t *out);

#endif
