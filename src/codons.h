/* (c) 2026 Giorgio Gonnella */

/*
   The genetic code: the translation of a coding sequence, used to find
   a conserved peptide in a reference (anchor-site.h) and to read the
   residues a read carries beside it (cleavage-site.h).

   The standard code, which influenza A uses. A codon with a character
   which is not A, C, G or T (an IUPAC code, an N) gives 'X', so that an
   ambiguous base is visible in the translation instead of being guessed;
   an incomplete codon at the end gives nothing, so that a translation is
   as long as the whole codons of the sequence.
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
   the first codon of <peptide>: the sequence is translated in its three
   frames and the first frame in which the peptide occurs gives the
   position. -1 if it occurs in none.
*/
long codons_find_peptide(const char *seq, size_t len, const char *peptide);

/* Unit tests */
void codons_test(void);

#endif
