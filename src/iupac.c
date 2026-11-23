/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <pthread.h>
#include <stdbool.h>
#include "iupac.h"

/*
   The pairs of characters which edlib must treat as equal: each IUPAC
   code with each of the bases it stands for, in both directions (a
   degenerate base can be in the primer, e.g. R, or in the read, e.g. N,
   which is also the padding character). Written out, so that the table
   is a constant: it is shared by the searchers of all threads (and by
   the scan of a sequence for a pattern, see pattern-scan.h).
*/
const EdlibEqualityPair iupac_equalities[] = {
  {'R', 'A'}, {'R', 'G'},
  {'Y', 'C'}, {'Y', 'T'},
  {'S', 'C'}, {'S', 'G'},
  {'W', 'A'}, {'W', 'T'},
  {'K', 'G'}, {'K', 'T'},
  {'M', 'A'}, {'M', 'C'},
  {'B', 'C'}, {'B', 'G'}, {'B', 'T'},
  {'D', 'A'}, {'D', 'G'}, {'D', 'T'},
  {'H', 'A'}, {'H', 'C'}, {'H', 'T'},
  {'V', 'A'}, {'V', 'C'}, {'V', 'G'},
  {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'},
  /* N (also the padding character) with the codes of the primer: a
     degenerate base of a primer overhanging the read is not an error */
  {'N', 'R'}, {'N', 'Y'}, {'N', 'S'}, {'N', 'W'}, {'N', 'K'}, {'N', 'M'},
  {'N', 'B'}, {'N', 'D'}, {'N', 'H'}, {'N', 'V'}, {'N', 'N'},
  {'U', 'T'},
};
const int n_iupac_equalities =
  (int)(sizeof(iupac_equalities) / sizeof(iupac_equalities[0]));

static iupac_equal_chars_t equal_chars[256];
static pthread_once_t equal_chars_once = PTHREAD_ONCE_INIT;

static void init_equal_chars(void) {
  for (int c = 0; c < 256; c++) {
    equal_chars[c].chars[0] = (unsigned char)c;
    equal_chars[c].n = 1;
  }
  for (int i = 0; i < n_iupac_equalities; i++) {
    unsigned char pair[2] = {(unsigned char)iupac_equalities[i].first,
                             (unsigned char)iupac_equalities[i].second};
    for (int side = 0; side < 2; side++) {
      iupac_equal_chars_t *a = equal_chars + pair[side];
      unsigned char b = pair[1 - side];
      bool known = false;
      for (int j = 0; j < a->n; j++)
        known = known || a->chars[j] == b;
      if (!known) {
        assert(a->n < IUPAC_MAX_EQUAL);
        a->chars[a->n++] = b;
      }
    }
  }
}

const iupac_equal_chars_t* iupac_equal_chars(void) {
  pthread_once(&equal_chars_once, init_equal_chars);
  return equal_chars;
}
