/* (c) 2026 Giorgio Gonnella */

/*
   The matches of the conserved termini of the primers of a scheme in a
   whole read, on both strands: where an amplicon can start or end
   inside a read (see chimera.h, which decides whether it does).

   On the strand of the read, a match is either the start of an amplicon
   (the terminus of a primer of the left end, or the reverse complement
   of that of a primer of the right end: the amplicon is then on the
   other strand) or the end of one (the other two forms).

   The termini are searched with the errors allowed in a primer, but at
   most TERMINUS_MATCHES_MAX_ERRORS: a terminus is short, and with more
   errors a read holds several matches by chance. The read is scanned
   once for each terminus (pattern-scan.h), so that all the matches are
   found; of matches which overlap, the one with the fewest errors is
   kept (then that of the first terminus of the scheme, the first in the
   read).

   A terminus is the inner part of a primer: a primer with a tail has
   the tail on the outer side of the match
   (terminus_matches_primer_bound()).
*/

#ifndef TERMINUS_MATCHES_H
#define TERMINUS_MATCHES_H

#include <stdbool.h>
#include <stddef.h>
#include "scheme.h"
#include "virus.h"

/* the most errors allowed in a match of a terminus */
#define TERMINUS_MATCHES_MAX_ERRORS 2

/* the greatest number of matches reported in a read */
#define TERMINUS_MATCHES_MAX 16

/* a match of a terminus in a read (0-based, inclusive) */
typedef struct {
  long from, to;
  long n_errors;
  size_t pattern;      /* the terminus matched (for this module) */
  bool start;          /* the start of an amplicon (else the end of one) */
} terminus_match_t;

typedef struct terminus_matches_t terminus_matches_t;

/*
   Prepare the search of the termini of <scheme> (it contains buffers:
   use one for each thread). <max_errors>: the errors allowed in a
   primer, a number of bases, or, if smaller than 1, a fraction of its
   length (in a terminus, at most TERMINUS_MATCHES_MAX_ERRORS).
*/
terminus_matches_t* terminus_matches_new(const scheme_t *scheme,
    double max_errors);

/*
   The matches of the termini in <seq> (length <len>), which do not
   overlap, in the order of the read, into <matches> (room for
   TERMINUS_MATCHES_MAX); their number is returned.
*/
size_t terminus_matches_find(terminus_matches_t *self, const char *seq,
    size_t len, terminus_match_t *matches);

/*
   Where the primer of the match <m> of <seq> (length <len>) begins (the
   start of an amplicon) or ends (the end of one): the primers with a
   tail which have the terminus are searched around the match, whole,
   with the errors allowed in a primer; without any (no tail, or not
   found), the bound of the match itself.
*/
long terminus_matches_primer_bound(terminus_matches_t *self,
    const char *seq, size_t len, const terminus_match_t *m);

/*
   Search only the termini of the primers of the types <types> (see
   virus.h): the termini of the other type are not in the reads of a
   sample of this one, and the short terminus of influenza A matches
   inside a sequence by chance.
*/
void terminus_matches_restrict_virus(terminus_matches_t *self,
    virus_set_t types);

/* Delete a prepared search (NULL is accepted) */
void terminus_matches_delete(terminus_matches_t *self);

/* Unit tests */
void terminus_matches_test(void);

#endif
