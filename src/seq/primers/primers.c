/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stddef.h>
#include <string.h>
/* the tables of the primers, written by the build from
   metadata/primers.yaml (scripts/primers-table.awk) */
#define PRIMERS_TABLE_DEFINE
#include "primers.h"
#include "iupac.h"
#include "sequtil.h"
#include "test.h"

const primer_t* primers_builtin(const char *name) {
  assert(name != NULL);
  for (size_t i = 0; i < PRIMERS_N_BUILTIN; i++) {
    const char *a = name, *b = PRIMERS_BUILTIN[i].name;
    while (*a != '\0' && *b != '\0' &&
           (*a == *b || (*a >= 'A' && *a <= 'Z' && *a + 32 == *b)))
      a++, b++;
    if (*a == '\0' && *b == '\0')
      return PRIMERS_BUILTIN + i;
  }
  return NULL;
}

const primer_t* primers_builtin_all(size_t *n) {
  if (n != NULL)
    *n = PRIMERS_N_BUILTIN;
  return PRIMERS_BUILTIN;
}

bool primers_binds_segment(const primer_t *p, int segment) {
  assert(p != NULL);
  if (p->segments == PRIMER_SEGMENTS_ALL)
    return true;
  if (segment < 1 || segment > SEGMENTS_N)
    return false;
  return (p->segments & PRIMER_SEGMENT(segment)) != 0;
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
  reverse_complement(seq, len, dest);
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
    test_assert(iupac_is_seq(all[i].seq, all[i].seq_len));
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
  /* the reverse complement, with degenerate codes */
  char rc[32];
  primers_reverse_complement("AGCRAAAGCAGG", 12, rc);
  test_eq_mem(rc, 12, "CCTGCTTTYGCT");
  primers_reverse_complement(PRIMERS_TERMINUS_13,
                             strlen(PRIMERS_TERMINUS_13), rc);
  test_eq_mem(rc, strlen(PRIMERS_TERMINUS_13), "CCTTGTTTCTACT");
}
