/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "codons.h"
#include "test.h"

/* the residues of the 64 codons, in the order of the bases below: the
   standard code, as the usual TCAG table */
static const char residues[] =
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

/* the index of a base in the order T, C, A, G; -1 for any other
   character */
static int base_index(char c) {
  switch (c) {
    case 'T': case 't': case 'U': case 'u': return 0;
    case 'C': case 'c': return 1;
    case 'A': case 'a': return 2;
    case 'G': case 'g': return 3;
    default: return -1;
  }
}

char codons_residue(const char *codon) {
  assert(codon != NULL);
  int a = base_index(codon[0]), b = base_index(codon[1]),
      c = base_index(codon[2]);
  if (a < 0 || b < 0 || c < 0)
    return 'X';
  return residues[16 * a + 4 * b + c];
}

size_t codons_translate(const char *seq, size_t len, char *out) {
  assert(seq != NULL || len == 0);
  assert(out != NULL);
  size_t n = len / 3;
  for (size_t i = 0; i < n; i++)
    out[i] = codons_residue(seq + 3 * i);
  out[n] = '\0';
  return n;
}

/* whether <protein> begins with <pattern>, in which '.' is any residue */
static bool matches(const char *protein, const char *pattern) {
  for (size_t i = 0; pattern[i] != '\0'; i++)
    if (protein[i] == '\0' ||
        (pattern[i] != '.' && pattern[i] != protein[i]))
      return false;
  return true;
}

long codons_find_peptide(const char *seq, size_t len, const char *pattern) {
  assert(pattern != NULL);
  if (seq == NULL || len < 3 || pattern[0] == '\0')
    return -1;
  char *protein;
  XMALLOC(protein, len / 3 + 1);
  long found = -1;
  for (size_t frame = 0; frame < 3 && found < 0; frame++) {
    if (len <= frame)
      break;
    size_t n = codons_translate(seq + frame, len - frame, protein);
    for (size_t i = 0; i < n && found < 0; i++)
      if (matches(protein + i, pattern))
        found = (long)frame + 3 * (long)i;
  }
  free(protein);
  return found;
}

/* --- unit tests --- */

void codons_test(void) {
  test_eq_char(codons_residue("ATG"), 'M');
  test_eq_char(codons_residue("atg"), 'M');
  test_eq_char(codons_residue("TAA"), '*');
  test_eq_char(codons_residue("TAG"), '*');
  test_eq_char(codons_residue("TGA"), '*');
  test_eq_char(codons_residue("GGN"), 'X');
  test_eq_char(codons_residue("TTT"), 'F');
  test_eq_char(codons_residue("GGG"), 'G');
  char out[16];
  /* GGA CTG TTT GGG GCT: the start of the fusion peptide */
  test_eq_size(codons_translate("GGACTGTTTGGGGCT", 15, out), 5);
  test_eq_str(out, "GLFGA");
  /* an incomplete codon at the end is not translated */
  test_eq_size(codons_translate("GGACTGTTTGGGGCTAT", 17, out), 5);
  test_eq_str(out, "GLFGA");
  test_eq_size(codons_translate("AT", 2, out), 0);
  test_eq_str(out, "");
  /* the peptide in each of the three frames */
  test_eq_long(codons_find_peptide("GGACTGTTTGGGGCT", 15, "GLFGA"), 0);
  test_eq_long(codons_find_peptide("AGGACTGTTTGGGGCT", 16, "GLFGA"), 1);
  test_eq_long(codons_find_peptide("AAGGACTGTTTGGGGCT", 17, "GLFGA"), 2);
  test_eq_long(codons_find_peptide("AAGGACTGTTTGGGGCT", 17, "GLFGAIAG"), -1);
  test_eq_long(codons_find_peptide(NULL, 0, "GLFGA"), -1);
  test_eq_long(codons_find_peptide("GGACTGTTTGGGGCT", 15, ""), -1);
  /* a '.' matches any residue: the start of HA2 is GLFGAIAG in most
     haemagglutinins, GIFGAIAG in some H3, GLLGAIAG and, in influenza B,
     GFFGAIAG */
  const char *ha2 = "GGACTGTTTGGGGCTATAGCAGGC";       /* GLFGAIAG */
  const char *h3 = "GGAATTTTTGGGGCTATAGCAGGC";        /* GIFGAIAG */
  const char *flub = "GGATTTTTTGGGGCTATAGCAGGC";      /* GFFGAIAG */
  test_eq_long(codons_find_peptide(ha2, 24, "G..GAIAG"), 0);
  test_eq_long(codons_find_peptide(h3, 24, "G..GAIAG"), 0);
  test_eq_long(codons_find_peptide(flub, 24, "G..GAIAG"), 0);
  test_eq_long(codons_find_peptide(h3, 24, "GLFGAIAG"), -1);
  /* the pattern must still fit entirely */
  test_eq_long(codons_find_peptide(ha2, 21, "G..GAIAG"), -1);
  /* a pattern of dots alone matches at the first codon */
  test_eq_long(codons_find_peptide(ha2, 24, "..."), 0);
}
