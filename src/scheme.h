/* (c) 2026 Giorgio Gonnella */

/*
   The primer scheme of a run: the primers which can be found at each end
   of an amplicon.

   An amplicon of the multisegment RT-PCR has, on its plus strand (the
   orientation of the mRNA), a primer binding the conserved 12 nt
   terminus at the left end and the reverse complement of a primer
   binding the conserved 13 nt terminus at the right end (see primers.h):

     [tail] AGCAAAAGCAGG ... segment ... CCTTGTTTCTACT [reverse tail]

   A scheme holds, for each end, the primers which may be found there:
   several, as a protocol can use variants of the same primer (e.g. the
   two forms of the Uni12 primer, which differ at the fourth position of
   the terminus). The primers of the right end are stored as they are
   found in a read, i.e. already reverse-complemented, so that both ends
   are searched in the same orientation.

   The conserved termini themselves are built in and cannot be changed:
   they are the binding sites on the segments. What a protocol changes is
   the tail added to them, which is why a scheme is given by the primers
   (tail and terminus together), not by the termini.
*/

#ifndef SCHEME_H
#define SCHEME_H

#include <stddef.h>
#include "logger.h"
#include "primers.h"

/* name of the default built-in scheme */
#define SCHEME_DEFAULT "uni12-mbt13"

typedef struct {
  char *name;
  /* the primers of the left end, as they are found in a read */
  primer_t *left;
  size_t n_left;
  /* the primers of the right end, reverse-complemented (as they are
     found in a read) */
  primer_t *right;
  size_t n_right;
  /* the sequences of the reverse-complemented primers of the right end */
  char *right_seqs;
} scheme_t;

/*
   Load a built-in scheme by name (case is ignored); the program is
   exited with an error message if there is none with that name.

   The built-in schemes:

   - "uni12-mbt13" (the default): the untailed degenerate Uni12 primer
     (uni12r) at the left end and the tailed mbtuni-13 at the right
     end. The two ends of
     a protocol need not use the same form of the primers; with an
     untailed Uni12 primer the left end of an amplicon is only the 12
     bases of the conserved terminus, which are short enough to anneal
     elsewhere in a segment (see amplicon.h on internal priming).
   - "mbt": the tailed primers of the multisegment RT-PCR
     of Zhou et al. 2009, with both forms of the Uni12 primer
     (mbtuni-12 and mbtuni-12.4r) at the left end and mbtuni-13 at the
     right end;
   - "uni": the conserved termini used as primers, without a tail
     (uni12 and uni12r at the left end, uni13 at the right end): for
     protocols which add no tail.
*/
scheme_t* scheme_new_builtin(const char *name, Logger *logger);

/* The names of the built-in schemes, NULL-terminated */
const char* const* scheme_builtin_names(void);

/* Delete a scheme (NULL is accepted) */
void scheme_delete(scheme_t *self);

/* Unit tests */
void scheme_test(void);

#endif
