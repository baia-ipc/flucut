/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stddef.h>
#include <string.h>
#include "primers.h"
#include "test.h"

/*
   The built-in primers (see primers.h for the sources of each sequence).

   The length of the oligonucleotide and that of its conserved terminus
   are computed by primers_builtin_init(), so that a sequence is written
   only once here.
*/
static primer_t builtin[] = {
  {"uni12", PRIMERS_TERMINUS_12, 0, 0, PRIMER_END_12},
  {"uni12r", PRIMERS_TERMINUS_12_R, 0, 0, PRIMER_END_12},
  {"uni13", PRIMERS_TERMINUS_13, 0, 0, PRIMER_END_13},
  {"mbtuni-12", PRIMERS_MBT_TAIL "CAAAAGCAGG", 0, 0, PRIMER_END_12},
  {"mbtuni-12.4r", PRIMERS_MBT_TAIL "CRAAAGCAGG", 0, 0, PRIMER_END_12},
  {"mbtuni-13", PRIMERS_MBT_TAIL "TAGAAACAAGG", 0, 0, PRIMER_END_13},
};
#define N_BUILTIN ((size_t)(sizeof(builtin) / sizeof(builtin[0])))

/* the conserved terminus bound by a primer of each end */
static const char* terminus_of(primer_end_t end) {
  return end == PRIMER_END_12 ? PRIMERS_TERMINUS_12 : PRIMERS_TERMINUS_13;
}

/*
   Compute the lengths of the built-in primers, once: the terminus is the
   suffix of the oligonucleotide of the length of the conserved terminus
   of its end (the tail and the terminus may overlap in the
   oligonucleotide, see primers.h)
*/
static void primers_builtin_init(void) {
  static bool done = false;
  if (done)
    return;
  for (size_t i = 0; i < N_BUILTIN; i++) {
    primer_t *p = builtin + i;
    p->seq_len = strlen(p->seq);
    p->terminus_len = strlen(terminus_of(p->end));
    assert(p->terminus_len <= p->seq_len);
  }
  done = true;
}

const primer_t* primers_builtin(const char *name) {
  assert(name != NULL);
  primers_builtin_init();
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
  primers_builtin_init();
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
