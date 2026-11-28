/* (c) 2026 Giorgio Gonnella */

/*
   Which characters of a primer, a probe or a read are equal: each IUPAC
   code is equal to each of the bases it stands for.
*/

#ifndef IUPAC_H
#define IUPAC_H

#include <stdbool.h>
#include <stddef.h>

/* A pair of characters that are equivalent for sequence matching. */
typedef struct { char first, second; } iupac_pair_t;

/* The IUPAC codes and bases they match, including symmetric pairs. */
extern const iupac_pair_t iupac_equalities[];
extern const size_t n_iupac_equalities;

/* The canonical bases represented by an IUPAC code, or NULL. */
const char* iupac_bases(char code);

/* Whether every character in the non-empty sequence is an IUPAC code. */
bool iupac_is_seq(const char *seq, size_t len);


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
