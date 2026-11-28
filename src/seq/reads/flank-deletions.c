/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "edlib.h"
#include "edlib-iupac.h"
#include "fasta.h"
#include "flank-deletions.h"
#include "test.h"

/* the bases of the reference beyond the window, for the deletions of the
   read: the reference is longer than the bases aligned to it */
#define FLANK_DELETIONS_SLACK 100

/* the greatest share of the bases aligned which may differ from the
   reference (as for the junction, see junction.c) */
#define FLANK_DELETIONS_MAX_DIFFERENCE 0.3

/* deletion runs separated by at most this many bases of the read are one
   deletion: an alignment by edit distance (no cost for opening a gap)
   can break a deletion into pieces, with a few bases between them */
#define FLANK_DELETIONS_MERGE 3

/*
   Align the <n> bases of <query> to <target> (<target_len> bases), from
   the first base of both on, the end of the target free: the number of
   bases of the query before each deletion of at least
   FLANK_DELETIONS_MIN bases, into <at> (at most <max>); their number. A
   deletion is a run of bases of the target not in the query (an
   insertion to the query in the terms of edlib), or several runs close
   together (FLANK_DELETIONS_MERGE); it counts only with bases of the
   query on both sides of it.
*/
static size_t deletions_in(const char *query, long n, const char *target,
    long target_len, long *at, size_t max) {
  int n_pairs;
  const EdlibEqualityPair *pairs = edlib_iupac_equalities(&n_pairs);
  EdlibAlignResult r = edlibAlign(query, (int)n, target, (int)target_len,
      edlibNewAlignConfig((int)(FLANK_DELETIONS_MAX_DIFFERENCE * (double)n),
                          EDLIB_MODE_SHW, EDLIB_TASK_PATH, pairs, n_pairs));
  size_t found = 0;
  if (r.status == EDLIB_STATUS_OK && r.editDistance >= 0) {
    /* the deletion being read: the bases of the query before it, its
       bases of the target, and the bases of the query read since its
       last run (-1: none being read) */
    long q = 0, start = -1, deleted = 0, since = 0;
    for (int i = 0; i < r.alignmentLength; i++) {
      if (r.alignment[i] == EDLIB_EDOP_DELETE) {
        if (start < 0 || since > FLANK_DELETIONS_MERGE) {
          if (start > 0 && deleted >= FLANK_DELETIONS_MIN && found < max)
            at[found++] = start;
          start = q;
          deleted = 0;
        }
        deleted++;
        since = 0;
        continue;
      }
      q++;
      since++;
    }
    /* (the last one, if bases of the query follow it; not one at the
       junction itself, whose bases belong to it) */
    if (start > 0 && deleted >= FLANK_DELETIONS_MIN && since > 0 &&
        found < max)
      at[found++] = start;
  }
  edlibFreeAlignResult(r);
  return found;
}

/* <n> bases of <from>, in the reverse order (not complemented) */
static void reversed_copy(const char *from, long n, char *to) {
  for (long i = 0; i < n; i++)
    to[i] = from[n - 1 - i];
}

size_t flank_deletions_find(const reference_set_t *refs, const char *seq,
    size_t len, size_t ref_index, const junction_t *j, long *cuts) {
  assert(refs != NULL);
  assert(j != NULL);
  assert(cuts != NULL);
  if (!j->found)
    return 0;
  const fasta_record_t *ref = reference_set_fasta(refs)->records + ref_index;
  const long n = (long)len, offset = j->offset;
  long at[FLANK_DELETIONS_MAX_CUTS];
  size_t found = 0;
  char *query, *target;
  XMALLOC(query, FLANK_DELETIONS_WINDOW);
  XMALLOC(target, FLANK_DELETIONS_WINDOW + FLANK_DELETIONS_SLACK);
  /* the left flank: its last bases and those of the reference before
     the junction, aligned backwards from it */
  long w = offset < FLANK_DELETIONS_WINDOW ? offset : FLANK_DELETIONS_WINDOW;
  long t_from = j->left_to - (w + FLANK_DELETIONS_SLACK);
  if (t_from < 0)
    t_from = 0;
  long t_len = j->left_to - t_from;
  if (w > 0 && t_len > 0 && (size_t)j->left_to <= ref->seq_len) {
    reversed_copy(seq + offset - w, w, query);
    reversed_copy(ref->seq + t_from, t_len, target);
    size_t k = deletions_in(query, w, target, t_len, at,
                            FLANK_DELETIONS_MAX_CUTS);
    /* (backwards: the nearest to the junction first) */
    for (size_t i = k; i > 0; i--)
      cuts[found++] = offset - at[i - 1];
  }
  /* the right flank: its first bases and those of the reference from
     the junction on, aligned forwards */
  w = n - offset < FLANK_DELETIONS_WINDOW ? n - offset :
      FLANK_DELETIONS_WINDOW;
  t_len = (long)ref->seq_len - j->right_from;
  if (t_len > w + FLANK_DELETIONS_SLACK)
    t_len = w + FLANK_DELETIONS_SLACK;
  if (w > 0 && t_len > 0 && found < FLANK_DELETIONS_MAX_CUTS) {
    size_t k = deletions_in(seq + offset, w, ref->seq + j->right_from,
                            t_len, at, FLANK_DELETIONS_MAX_CUTS - found);
    for (size_t i = 0; i < k; i++)
      cuts[found++] = offset + at[i];
  }
  free(query);
  free(target);
  return found;
}

/* --- unit tests --- */

void flank_deletions_test(void) {
  enum { REF_LEN = 1000, PROBE = 32 };
  unsigned long state = 7;
  char ref[REF_LEN + 1], read[REF_LEN + 1];
  test_random_dna_seq(ref, REF_LEN, &state);
  fasta_record_t rec = {"r0", "", ref, REF_LEN};
  fasta_t fa = {&rec, 1};
  reference_index_t *ix = reference_index_new(&fa);
  reference_set_t *rs = reference_set_new(ix, PROBE, 0.2);
  long cuts[FLANK_DELETIONS_MAX_CUTS];

  /* a junction (the bases 300 to 699 missing), and nothing else: no
     cut */
  memcpy(read, ref, 300);
  memcpy(read + 300, ref + 700, 300);
  junction_t j = {.found = true, .left_to = 300, .right_from = 700,
                  .offset = 300};
  test_eq_size(flank_deletions_find(rs, read, 600, 0, &j, cuts), 0);
  /* not located: no cut */
  j.found = false;
  test_eq_size(flank_deletions_find(rs, read, 600, 0, &j, cuts), 0);
  j.found = true;

  /* 20 bases missing 100 bases before the junction (the bases 180 to
     199), and 8 bases missing 60 bases after it (the bases 760 to 767):
     a cut at each, at the first base after it in the template */
  size_t n = 0;
  memcpy(read, ref, 180);
  n = 180;
  memcpy(read + n, ref + 200, 100);   /* to the junction */
  n += 100;
  long offset = (long)n;
  memcpy(read + n, ref + 700, 60);
  n += 60;
  memcpy(read + n, ref + 768, 232);
  n += 232;
  j.offset = offset;
  size_t k = flank_deletions_find(rs, read, n, 0, &j, cuts);
  test_eq_size(k, 2);
  test_eq_long(cuts[0], 180);
  /* (the last base before the deletion is also the first deleted: the
     deletion can be put one base later) */
  test_assert(cuts[1] == offset + 60 || cuts[1] == offset + 61);

  /* a deletion shorter than FLANK_DELETIONS_MIN, or farther from the
     junction than FLANK_DELETIONS_WINDOW: no cut */
  n = 0;
  memcpy(read, ref, 20);
  n = 20;
  memcpy(read + n, ref + 50, 250);    /* 30 bases missing, 280 before */
  n += 250;
  offset = (long)n;
  memcpy(read + n, ref + 700, 100);
  n += 100;
  memcpy(read + n, ref + 803, 197);   /* 3 bases missing */
  n += 197;
  j.offset = offset;
  test_eq_size(flank_deletions_find(rs, read, n, 0, &j, cuts), 0);

  reference_set_delete(rs);
  reference_index_delete(ix);
}
