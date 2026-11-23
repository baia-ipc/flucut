/* (c) 2026 Giorgio Gonnella */

/*
   Which characters of a primer, a probe or a read are equal: each IUPAC
   code is equal to each of the bases it stands for.
*/

#ifndef IUPAC_H
#define IUPAC_H

#include "edlib.h"

/* the pairs of characters edlib must treat as equal: the IUPAC codes and
   the bases they stand for (shared by all edlib alignments of flucut) */
extern const EdlibEqualityPair iupac_equalities[];
extern const int n_iupac_equalities;

/* the most characters equal to one character */
#define IUPAC_MAX_EQUAL 24

/* the <n> characters equal to a character (itself first) */
typedef struct {
  unsigned char n;
  unsigned char chars[IUPAC_MAX_EQUAL];
} iupac_equal_chars_t;

/* For each of the 256 characters, the characters equal to it, from
   iupac_equalities (computed at the first call; safe with threads) */
const iupac_equal_chars_t* iupac_equal_chars(void);

#endif
