/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stddef.h>
#include <string.h>
#include "primers.h"
#include "test.h"

/* length of a string literal */
#define LITLEN(S) (sizeof(S) - 1)

/*
   A built-in primer: its sequence <SEQ> is the oligonucleotide, whose
   last bases are the conserved terminus of the end <END> (the lengths
   are computed by the compiler, so that each sequence is written once)
*/
#define BUILTIN(NAME, SEQ, TERMINUS, END) \
  {NAME, SEQ, LITLEN(SEQ), LITLEN(TERMINUS), END}

/* The built-in primers (see primers.h for the source of each sequence) */
static const primer_t builtin[] = {
  BUILTIN("uni12", PRIMERS_TERMINUS_12, PRIMERS_TERMINUS_12,
          PRIMER_END_12),
  BUILTIN("uni12r", PRIMERS_TERMINUS_12_R, PRIMERS_TERMINUS_12_R,
          PRIMER_END_12),
  BUILTIN("uni13", PRIMERS_TERMINUS_13, PRIMERS_TERMINUS_13,
          PRIMER_END_13),
  BUILTIN("mbtuni-12", PRIMERS_MBT_TAIL "CAAAAGCAGG",
          PRIMERS_TERMINUS_12, PRIMER_END_12),
  BUILTIN("mbtuni-12.4r", PRIMERS_MBT_TAIL "CRAAAGCAGG",
          PRIMERS_TERMINUS_12_R, PRIMER_END_12),
  BUILTIN("mbtuni-13", PRIMERS_MBT_TAIL "TAGAAACAAGG",
          PRIMERS_TERMINUS_13, PRIMER_END_13),
};
#define N_BUILTIN ((size_t)(sizeof(builtin) / sizeof(builtin[0])))

const primer_t* primers_builtin(const char *name) {
  assert(name != NULL);
  for (size_t i = 0; i < N_BUILTIN; i++) {
    const char *a = name, *b = builtin[i].name;
    while (*a != '\0' && *b != '\0' &&
           (*a == *b || (*a >= 'A' && *a <= 'Z' && *a + 32 == *b)))
      a++, b++;
    if (*a == '\0' && *b == '\0')
      return builtin + i;
  }
  return NULL;
}

const primer_t* primers_builtin_all(size_t *n) {
  if (n != NULL)
    *n = N_BUILTIN;
  return builtin;
}

/* the IUPAC codes and the bases they stand for */
static const struct { char code; const char *bases; } iupac[] = {
  {'A', "A"}, {'C', "C"}, {'G', "G"}, {'T', "T"}, {'U', "T"},
  {'R', "AG"}, {'Y', "CT"}, {'S', "CG"}, {'W', "AT"}, {'K', "GT"},
  {'M', "AC"}, {'B', "CGT"}, {'D', "AGT"}, {'H', "ACT"}, {'V', "ACG"},
  {'N', "ACGT"},
};
#define N_IUPAC ((size_t)(sizeof(iupac) / sizeof(iupac[0])))

const char* primers_iupac_bases(char code) {
  if (code >= 'a' && code <= 'z')
    code = (char)(code - 32);
  for (size_t i = 0; i < N_IUPAC; i++)
    if (iupac[i].code == code)
      return iupac[i].bases;
  return NULL;
}

bool primers_is_iupac_seq(const char *seq, size_t len) {
  assert(seq != NULL);
  if (len == 0)
    return false;
  for (size_t i = 0; i < len; i++)
    if (primers_iupac_bases(seq[i]) == NULL)
      return false;
  return true;
}

/* the complement of a IUPAC code (the code of the complemented bases) */
static char complement(char code) {
  static const char *from = "ACGTURYSWKMBDHVN";
  static const char *to   = "TGCAAYRSWMKVHDBN";
  if (code >= 'a' && code <= 'z')
    code = (char)(code - 32);
  const char *p = strchr(from, code);
  assert(p != NULL);
  return to[p - from];
}

void primers_reverse_complement(const char *seq, size_t len, char *dest) {
  assert(seq != NULL);
  assert(dest != NULL);
  for (size_t i = 0; i < len; i++)
    dest[i] = complement(seq[len - 1 - i]);
}

/* the terminus of <p> is the suffix of its sequence, and the part before
   it is its tail */
static void check_builtin(const char *name, const char *seq,
    const char *terminus, const char *tail) {
  const primer_t *p = primers_builtin(name);
  test_assert(p != NULL);
  test_eq_str(p->seq, seq);
  test_eq_size(p->seq_len, strlen(seq));
  test_eq_size(p->terminus_len, strlen(terminus));
  test_eq_mem(p->seq + p->seq_len - p->terminus_len, p->terminus_len,
              terminus);
  test_eq_mem(p->seq, p->seq_len - p->terminus_len, tail);
}

void primers_test(void) {
  /* the sequences of the publications (see primers.h) */
  check_builtin("uni12", "AGCAAAAGCAGG", "AGCAAAAGCAGG", "");
  check_builtin("uni13", "AGTAGAAACAAGG", "AGTAGAAACAAGG", "");
  check_builtin("uni12r", "AGCRAAAGCAGG", "AGCRAAAGCAGG", "");
  check_builtin("mbtuni-12", "ACGCGTGATCAGCAAAAGCAGG", "AGCAAAAGCAGG",
                "ACGCGTGATC");
  check_builtin("mbtuni-13", "ACGCGTGATCAGTAGAAACAAGG", "AGTAGAAACAAGG",
                "ACGCGTGATC");
  check_builtin("mbtuni-12.4r", "ACGCGTGATCAGCRAAAGCAGG", "AGCRAAAGCAGG",
                "ACGCGTGATC");
  /* the published tail is the start of both tailed primers: its last two
     bases are also the first two of the conserved terminus */
  const primer_t *p = primers_builtin("mbtuni-12");
  test_eq_mem(p->seq, strlen(PRIMERS_MBT_TAIL), PRIMERS_MBT_TAIL);
  test_eq_size(strlen(PRIMERS_MBT_TAIL) + p->terminus_len - p->seq_len, 2);
  /* the names are matched without case */
  test_assert(primers_builtin("MBTuni-12") == p);
  test_assert(primers_builtin("mbtuni-14") == NULL);
  /* the primers of each end bind the same terminus */
  size_t n;
  const primer_t *all = primers_builtin_all(&n);
  test_eq_size(n, 6);
  for (size_t i = 0; i < n; i++)
    test_assert(primers_is_iupac_seq(all[i].seq, all[i].seq_len));
  /* IUPAC codes */
  test_eq_str(primers_iupac_bases('R'), "AG");
  test_eq_str(primers_iupac_bases('r'), "AG");
  test_eq_str(primers_iupac_bases('A'), "A");
  test_assert(primers_iupac_bases('X') == NULL);
  test_assert(!primers_is_iupac_seq("ACGX", 4));
  test_assert(!primers_is_iupac_seq("", 0));
  /* the reverse complement, with degenerate codes */
  char rc[32];
  primers_reverse_complement("AGCRAAAGCAGG", 12, rc);
  test_eq_mem(rc, 12, "CCTGCTTTYGCT");
  primers_reverse_complement(PRIMERS_TERMINUS_13,
                             strlen(PRIMERS_TERMINUS_13), rc);
  test_eq_mem(rc, strlen(PRIMERS_TERMINUS_13), "CCTTGTTTCTACT");
}
