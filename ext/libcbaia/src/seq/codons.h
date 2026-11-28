/* (c) 2026 Giorgio Gonnella */

/*
   Translation with the standard genetic code. A codon with a character
   which is not A, C, G or T (for example an ambiguous base) gives 'X';
   an incomplete codon at the end is not translated.
*/

#ifndef CODONS_H
#define CODONS_H

#include <stddef.h>

/* The residue of the three bases at <codon>, upper or lower case: a
   letter, '*' for a stop codon, 'X' for a codon with another character */
char codons_residue(const char *codon);

/*
   Translate the <len> bases of <seq> into <out>, which must have room
   for len / 3 + 1 characters; the result is terminated. Returns the
   number of residues written (len / 3).
*/
size_t codons_translate(const char *seq, size_t len, char *out);

/*
   The 0-based position in <seq>, of <len> bases, of the first base of
   the first codon which matches <pattern>: the sequence is translated
   in its three frames and the first frame in which the pattern occurs
   gives the position. -1 if it occurs in none.

   A '.' in the pattern matches any residue, so that a peptide which is
   conserved at some of its positions only can still be looked for.
*/
long codons_find_peptide(const char *seq, size_t len, const char *pattern);

/* Unit tests */
void codons_test(void);

#endif
