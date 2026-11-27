/* (c) 2026 Giorgio Gonnella */

/*
   The matches of the primers of a scheme inside a read: between the
   primers of its two ends, in the part of the read which is kept (see
   amplicon.h: internal priming and chimeras).

   All the primers of both ends of the scheme are searched, whole, in
   the orientation in which the layout of the read was found, with the
   errors allowed in each. The read is scanned once for each primer
   (pattern-scan.h), so that all the matches are found; of matches which
   overlap, the one with the fewest errors is kept (then that of a
   primer of the left end, the longer one, the first primer of the
   scheme, the first in the read). The matches are not trimmed.

   A prepared search contains buffers: use one for each thread.
*/

#ifndef INTERNAL_MATCHES_H
#define INTERNAL_MATCHES_H

#include <stddef.h>
#include "amplicon.h"
#include "scheme.h"
#include "virus.h"

typedef struct internal_matches_t internal_matches_t;

/*
   Prepare the search of the primers of <scheme>, which must exist as
   long as the result is used; <left_max_errors>, <right_max_errors>:
   the errors allowed in each primer of each end.
*/
internal_matches_t* internal_matches_new(const scheme_t *scheme,
    const long *left_max_errors, const long *right_max_errors);

/*
   The matches in the bases <from> to <to> - 1 of <seq>, in the order of
   the read, into <matches> (room for AMPLICON_MAX_INTERNAL, the
   positions are those in <seq>); their number is returned.
*/
size_t internal_matches_find(internal_matches_t *self, const char *seq,
    long from, long to, amplicon_match_t *matches);

/*
   Search only the primers of the types <types> (see virus.h): the
   primers of the other type were not in the reaction of a sample of
   this one, and the short primers of influenza A match inside a
   sequence by chance, which would be reported as sites of internal
   priming.
*/
void internal_matches_restrict_virus(internal_matches_t *self,
    virus_set_t types);

/* Delete a prepared search (NULL is accepted) */
void internal_matches_delete(internal_matches_t *self);

/* Unit tests */
void internal_matches_test(void);

#endif
