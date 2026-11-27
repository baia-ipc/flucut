/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
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

long codons_find_peptide(const char *seq, size_t len, const char *peptide) {
  assert(peptide != NULL);
  if (seq == NULL || len < 3 || peptide[0] == '\0')
    return -1;
  char *protein;
  XMALLOC(protein, len / 3 + 1);
  long found = -1;
  for (size_t frame = 0; frame < 3 && found < 0; frame++) {
    if (len <= frame)
      break;
    codons_translate(seq + frame, len - frame, protein);
    const char *hit = strstr(protein, peptide);
    if (hit != NULL)
      found = (long)frame + 3 * (long)(hit - protein);
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
}
