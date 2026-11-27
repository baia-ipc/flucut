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
   are computed by the compiler, so that each sequence is written once);
   <VIRUS> is the type of influenza it was designed for and <SEGMENTS>
   the segments whose terminus it binds (see primers.h)
*/
#define BUILTIN(NAME, SEQ, TERMINUS, END, VIRUS, SEGMENTS) \
  {NAME, SEQ, LITLEN(SEQ), LITLEN(SEQ) - LITLEN(TERMINUS), \
   LITLEN(TERMINUS), END, VIRUS, SEGMENTS}

/* a universal primer of influenza A: it binds every segment */
#define BUILTIN_A(NAME, SEQ, TERMINUS, END) \
  BUILTIN(NAME, SEQ, TERMINUS, END, VIRUS_A, PRIMER_SEGMENTS_ALL)

/*
   A primer of influenza B: <TAIL> is its 5' tail and <TERMINUS> the
   rest, which is the sequence of the segments <SEGMENTS> at the end
   <END> (see primers.h)
*/
#define BUILTIN_B(NAME, TAIL, TERMINUS, END, SEGMENTS) \
  BUILTIN(NAME, TAIL TERMINUS, TERMINUS, END, VIRUS_B, SEGMENTS)

/* the left and right ends of a segment of influenza B, as they are in
   the primers: the conserved terminus and the bases next to it */
#define IBV_LEFT(NEXT) PRIMERS_TERMINUS_B_9 NEXT
#define IBV_RIGHT(NEXT) PRIMERS_TERMINUS_B_10 NEXT
/* the form of the 10 nt terminus of the HA, the NA and the NS */
#define IBV_RIGHT_T(NEXT) PRIMERS_TERMINUS_B_10_T NEXT

/* The built-in primers (see primers.h for the source of each sequence) */
static const primer_t builtin[] = {
  /* influenza A (Uni12/Uni13) */
  BUILTIN_A("uni12", PRIMERS_TERMINUS_12, PRIMERS_TERMINUS_12,
            PRIMER_END_LEFT),
  BUILTIN_A("uni12r", PRIMERS_TERMINUS_12_R, PRIMERS_TERMINUS_12_R,
            PRIMER_END_LEFT),
  BUILTIN_A("uni13", PRIMERS_TERMINUS_13, PRIMERS_TERMINUS_13,
            PRIMER_END_RIGHT),
  BUILTIN_A("mbtuni-12", PRIMERS_MBT_TAIL "CAAAAGCAGG",
            PRIMERS_TERMINUS_12, PRIMER_END_LEFT),
  BUILTIN_A("mbtuni-12.4r", PRIMERS_MBT_TAIL "CRAAAGCAGG",
            PRIMERS_TERMINUS_12_R, PRIMER_END_LEFT),
  BUILTIN_A("mbtuni-13", PRIMERS_MBT_TAIL "TAGAAACAAGG",
            PRIMERS_TERMINUS_13, PRIMER_END_RIGHT),
  /* influenza B (the IBV-GA2 cocktail), each with the segments it
     binds (see primers.h) */
  BUILTIN_B("b-pbs-unif", PRIMERS_IBV_TAIL_LEFT, IBV_LEFT("GGAGC"),
            PRIMER_END_LEFT,
            PRIMER_SEGMENT(SEGMENT_PB1) | PRIMER_SEGMENT(SEGMENT_PB2)),
  BUILTIN_B("b-pbs-unir", PRIMERS_IBV_TAIL_RIGHT, IBV_RIGHT("CGAGC"),
            PRIMER_END_RIGHT,
            PRIMER_SEGMENT(SEGMENT_PB1) | PRIMER_SEGMENT(SEGMENT_PB2)),
  BUILTIN_B("b-pa-unif", PRIMERS_IBV_TAIL_LEFT, IBV_LEFT("GGTGC"),
            PRIMER_END_LEFT, PRIMER_SEGMENT(SEGMENT_PA)),
  BUILTIN_B("b-pa-unir", PRIMERS_IBV_TAIL_RIGHT, IBV_RIGHT("CGTGC"),
            PRIMER_END_RIGHT, PRIMER_SEGMENT(SEGMENT_PA)),
  BUILTIN_B("b-hana-unif", PRIMERS_IBV_TAIL_LEFT, IBV_LEFT("AGAGC"),
            PRIMER_END_LEFT,
            PRIMER_SEGMENT(SEGMENT_HA) | PRIMER_SEGMENT(SEGMENT_NA)),
  BUILTIN_B("b-hana-unir", PRIMERS_IBV_TAIL_RIGHT, IBV_RIGHT_T("AGAGC"),
            PRIMER_END_RIGHT,
            PRIMER_SEGMENT(SEGMENT_HA) | PRIMER_SEGMENT(SEGMENT_NA)),
  BUILTIN_B("b-np-unif", PRIMERS_IBV_TAIL_LEFT, IBV_LEFT("ACAGC"),
            PRIMER_END_LEFT, PRIMER_SEGMENT(SEGMENT_NP)),
  BUILTIN_B("b-np-unir", PRIMERS_IBV_TAIL_RIGHT, IBV_RIGHT("ACAGC"),
            PRIMER_END_RIGHT, PRIMER_SEGMENT(SEGMENT_NP)),
  BUILTIN_B("b-m-uni3f", PRIMERS_IBV_TAIL_LEFT, IBV_LEFT("ACGCACTT"),
            PRIMER_END_LEFT, PRIMER_SEGMENT(SEGMENT_M)),
  BUILTIN_B("b-mg-uni3f", PRIMERS_IBV_TAIL_LEFT, IBV_LEFT("AGGCACTT"),
            PRIMER_END_LEFT, PRIMER_SEGMENT(SEGMENT_M)),
  BUILTIN_B("b-m-uni3r", PRIMERS_IBV_TAIL_RIGHT, IBV_RIGHT("ACGCACTT"),
            PRIMER_END_RIGHT, PRIMER_SEGMENT(SEGMENT_M)),
  BUILTIN_B("b-ns-uni3f", PRIMERS_IBV_TAIL_LEFT, IBV_LEFT("AGAGGATT"),
            PRIMER_END_LEFT, PRIMER_SEGMENT(SEGMENT_NS)),
  BUILTIN_B("b-ns-uni3r", PRIMERS_IBV_TAIL_RIGHT, IBV_RIGHT_T("AGAGGATT"),
            PRIMER_END_RIGHT, PRIMER_SEGMENT(SEGMENT_NS)),
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

bool primers_binds_segment(const primer_t *p, int segment) {
  assert(p != NULL);
  if (p->segments == PRIMER_SEGMENTS_ALL)
    return true;
  if (segment < 1 || segment > SEGMENTS_N)
    return false;
  return (p->segments & PRIMER_SEGMENT(segment)) != 0;
}

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

void primers_reverse(const primer_t *p, char *seq_buf, primer_t *dest) {
  assert(p != NULL);
  assert(seq_buf != NULL);
  assert(dest != NULL);
  primers_reverse_complement(p->seq, p->seq_len, seq_buf);
  dest->name = p->name;
  dest->seq = seq_buf;
  dest->seq_len = p->seq_len;
  dest->terminus_from = p->seq_len - (p->terminus_from + p->terminus_len);
  dest->terminus_len = p->terminus_len;
  dest->end = p->end;
  dest->virus = p->virus;
  dest->segments = p->segments;
}

void primers_reverse_complement(const char *seq, size_t len, char *dest) {
  assert(seq != NULL);
  assert(dest != NULL);
  for (size_t i = 0; i < len; i++)
    dest[i] = complement(seq[len - 1 - i]);
}

void primers_test(void) {
  /* the published tail is the start of both tailed primers: its last two
     bases are also the first two of the conserved terminus */
  const primer_t *p = primers_builtin("mbtuni-12");
  test_eq_mem(p->seq, strlen(PRIMERS_MBT_TAIL), PRIMERS_MBT_TAIL);
  test_eq_size(strlen(PRIMERS_MBT_TAIL) + p->terminus_len - p->seq_len, 2);
  /* the names are matched without case */
  test_assert(primers_builtin("MBTuni-12") == p);
  test_assert(primers_builtin("mbtuni-14") == NULL);
  /* the built-in sequences contain only IUPAC codes, and the terminus
     of each lies inside it */
  size_t n;
  const primer_t *all = primers_builtin_all(&n);
  for (size_t i = 0; i < n; i++) {
    test_assert(primers_is_iupac_seq(all[i].seq, all[i].seq_len));
    test_assert(all[i].terminus_from + all[i].terminus_len <=
                all[i].seq_len);
  }
  /* the universal primers of influenza A bind every segment */
  test_assert(primers_binds_segment(primers_builtin("uni12r"), 4));
  test_assert(primers_binds_segment(primers_builtin("mbtuni-13"), 1));
  /* a primer of influenza B binds the segments it was designed for and
     no others, and an unknown segment is not one of them */
  const primer_t *b = primers_builtin("b-hana-unir");
  test_eq_int((int)b->virus, (int)VIRUS_B);
  test_assert(primers_binds_segment(b, SEGMENT_HA));
  test_assert(primers_binds_segment(b, SEGMENT_NA));
  test_assert(!primers_binds_segment(b, SEGMENT_PB1));
  test_assert(!primers_binds_segment(b, 0));
  test_assert(primers_binds_segment(primers_builtin("b-pbs-unif"),
                                    SEGMENT_PB2));
  test_assert(!primers_binds_segment(primers_builtin("b-m-uni3f"),
                                     SEGMENT_NS));
  /* each primer of influenza B is its tail and then the sequence of
     the segments it binds: the conserved terminus of its end and the
     bases next to it (5 or 8) */
  for (size_t i = 0; i < n; i++) {
    if (all[i].virus != VIRUS_B)
      continue;
    const char *tail = all[i].end == PRIMER_END_LEFT ?
        PRIMERS_IBV_TAIL_LEFT : PRIMERS_IBV_TAIL_RIGHT;
    test_eq_size(all[i].terminus_from, strlen(tail));
    test_eq_mem(all[i].seq, strlen(tail), tail);
    const char *terminus = all[i].seq + all[i].terminus_from;
    size_t universal = all[i].end == PRIMER_END_LEFT ?
        strlen(PRIMERS_TERMINUS_B_9) : strlen(PRIMERS_TERMINUS_B_10);
    test_assert(all[i].terminus_len == universal + 5 ||
                all[i].terminus_len == universal + 8);
    if (all[i].end == PRIMER_END_LEFT)
      test_eq_mem(terminus, universal, PRIMERS_TERMINUS_B_9);
    else
      test_assert(memcmp(terminus, PRIMERS_TERMINUS_B_10, universal) == 0
          || memcmp(terminus, PRIMERS_TERMINUS_B_10_T, universal) == 0);
    /* the segments it binds are those of its name, never none */
    test_assert(all[i].segments != PRIMER_SEGMENTS_ALL);
  }
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
