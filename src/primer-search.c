/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "edlib.h"
#include "alloc.h"
#include "numutils.h"
#include "primer-search.h"
#include "primers.h"
#include "test.h"

/* the character the read is padded with, to allow a primer to overhang
   it; it matches any base (it is the IUPAC code for "any base") */
#define PAD_CHAR 'N'

/*
   The pairs of characters which edlib must treat as equal: each IUPAC
   code with each of the bases it stands for, in both directions (a
   degenerate base can be in the primer, e.g. R, or in the read, e.g. N,
   which is also the padding character). Written out, so that the table
   is a constant: it is shared by the searchers of all threads.
*/
static const EdlibEqualityPair iupac_equalities[] = {
  {'R', 'A'}, {'R', 'G'},
  {'Y', 'C'}, {'Y', 'T'},
  {'S', 'C'}, {'S', 'G'},
  {'W', 'A'}, {'W', 'T'},
  {'K', 'G'}, {'K', 'T'},
  {'M', 'A'}, {'M', 'C'},
  {'B', 'C'}, {'B', 'G'}, {'B', 'T'},
  {'D', 'A'}, {'D', 'G'}, {'D', 'T'},
  {'H', 'A'}, {'H', 'C'}, {'H', 'T'},
  {'V', 'A'}, {'V', 'C'}, {'V', 'G'},
  {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'},
  {'U', 'T'},
};
#define N_IUPAC_EQUALITIES \
  ((int)(sizeof(iupac_equalities) / sizeof(iupac_equalities[0])))

struct primer_search_t {
  char *padded;
  size_t padded_cap;
};

primer_search_t* primer_search_new(void) {
  primer_search_t *self;
  XMALLOC(self, sizeof(primer_search_t));
  self->padded = NULL;
  self->padded_cap = 0;
  return self;
}

void primer_search_delete(primer_search_t *self) {
  if (self == NULL)
    return;
  free(self->padded);
  free(self);
}

/* edit distance between the primer and target[start..end] */
static long distance(const char *primer, size_t primer_len,
    const char *target, long start, long end) {
  EdlibAlignResult result = edlibAlign(primer, (int)primer_len,
      target + start, (int)(end - start + 1),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE,
                          iupac_equalities, N_IUPAC_EQUALITIES));
  assert(result.status == EDLIB_STATUS_OK);
  long d = result.editDistance;
  edlibFreeAlignResult(result);
  return d;
}

/*
   Positions in <target> of the bases of the primer at positions
   <pos_a> <= <pos_b>, for the alignment of the whole primer to
   target[start..end] (the alignment path is computed only for a primer
   which was found).

   A primer base aligned to a base of the target gives its position; a
   primer base which is deleted in the target gives the position of the
   target base which follows the deletion.

   The operations of the alignment path of edlib consume, in the
   orientation used here (the primer is the query, the read the target):
   EDLIB_EDOP_MATCH and EDLIB_EDOP_MISMATCH a base of each,
   EDLIB_EDOP_INSERT a base of the primer only, EDLIB_EDOP_DELETE a base
   of the read only (the names are those of edlib, which calls an
   insertion what is seen as a deletion from the read; the unit tests
   below check the resulting positions on reads with an insertion and
   with a deletion).
*/
static void map_primer_positions(const char *primer, size_t primer_len,
    const char *target, long start, long end, size_t pos_a, size_t pos_b,
    long *out_a, long *out_b) {
  assert(pos_a <= pos_b);
  assert(pos_b <= primer_len);
  *out_a = start;
  *out_b = end + 1;
  if (pos_a == 0 && pos_b == 0)
    return;
  EdlibAlignResult result = edlibAlign(primer, (int)primer_len,
      target + start, (int)(end - start + 1),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH,
                          iupac_equalities, N_IUPAC_EQUALITIES));
  assert(result.status == EDLIB_STATUS_OK);
  long qpos = 0, tpos = start;
  bool got_a = pos_a == 0, got_b = false;
  for (int i = 0; i < result.alignmentLength && !got_b; i++) {
    /* a base of the read inserted at a sought position belongs to the
       part of the primer before it (it is consumed before the position
       is reached): an inserted base at the end of the tail is tail, and
       is thus removed, rather than kept as a base of the terminus */
    if (result.alignment[i] == EDLIB_EDOP_DELETE) {
      tpos++;
      continue;
    }
    if (!got_a && (size_t)qpos == pos_a) {
      *out_a = tpos;
      got_a = true;
    }
    if ((size_t)qpos == pos_b) {
      *out_b = tpos;
      got_b = true;
      break;
    }
    switch (result.alignment[i]) {
      case EDLIB_EDOP_MATCH:
      case EDLIB_EDOP_MISMATCH:
        qpos++;
        tpos++;
        break;
      case EDLIB_EDOP_INSERT:    /* a base of the primer, not of the read */
        qpos++;
        break;
      default:
        assert(false);
    }
  }
  /* a position at the end of the primer is reached after the last
     operation */
  if (!got_a && (size_t)qpos == pos_a)
    *out_a = tpos;
  if (!got_b && (size_t)qpos == pos_b)
    *out_b = tpos;
  edlibFreeAlignResult(result);
}

/*
   Search the primer in <target> (the read, of length <read_len>, padded
   by <pad_left> characters on the left and <target_len> - <read_len> -
   <pad_left> on the right); returns false if not found; otherwise the
   positions in the read (clipped to the read), whether the primer
   overhangs the read, the number of errors and the number of primer
   positions inside the read (in <in_read>, unless NULL)
*/
static bool search(const char *primer, size_t primer_len,
    size_t terminus_from, size_t terminus_len, const char *target,
    size_t target_len, long pad_left, long read_len, long max_errors,
    primer_hit_t *hit, long *in_read) {
  if (target_len == 0)
    return false;
  EdlibAlignResult result = edlibAlign(primer, (int)primer_len,
      target, (int)target_len,
      edlibNewAlignConfig((int)max_errors, EDLIB_MODE_HW, EDLIB_TASK_LOC,
                          iupac_equalities, N_IUPAC_EQUALITIES));
  assert(result.status == EDLIB_STATUS_OK);
  long n_errors = result.editDistance;
  if (n_errors < 0) {
    edlibFreeAlignResult(result);
    return false;
  }
  /* among equally good locations, prefer the (first) longest one */
  long start = result.startLocations[0], end = result.endLocations[0];
  for (int i = 1; i < result.numLocations; i++)
    if (result.endLocations[i] - result.startLocations[i] > end - start) {
      start = result.startLocations[i];
      end = result.endLocations[i];
    }
  edlibFreeAlignResult(result);
  /* extend the aligned region while the edit distance does not increase */
  while (start > 0 &&
         distance(primer, primer_len, target, start - 1, end) <= n_errors)
    start--;
  while (end < (long)target_len - 1 &&
         distance(primer, primer_len, target, start, end + 1) <= n_errors)
    end++;
  /* the bounds of the conserved terminus, before clipping the alignment
     to the read */
  long terminus_start, terminus_end;
  map_primer_positions(primer, primer_len, target, start, end,
                       terminus_from, terminus_from + terminus_len,
                       &terminus_start, &terminus_end);
  hit->partial = start < pad_left || end >= pad_left + read_len;
  start -= pad_left;
  end -= pad_left;
  long clipped_start = start > 0 ? start : 0;
  long clipped_end = end < read_len - 1 ? end : read_len - 1;
  if (in_read != NULL)
    *in_read = clipped_end - clipped_start + 1;
  hit->start = clipped_start;
  hit->end = clipped_end;
  hit->n_errors = n_errors;
  hit->terminus_start = terminus_start - pad_left;
  hit->terminus_end = terminus_end - pad_left;
  return true;
}

static void not_found(primer_hit_t *hit) {
  hit->start = -1;
  hit->end = -1;
  hit->partial = false;
  hit->n_errors = -1;
  hit->terminus_start = -1;
  hit->terminus_end = -1;
}

/* the search of a sequence whose last <terminus_len> bases are the
   conserved terminus (0: none) */
static bool run(primer_search_t *self, const char *seq, size_t seq_len,
    size_t terminus_from, size_t terminus_len, const char *read,
    size_t read_len, long max_errors, bool partial_left,
    bool partial_right, primer_hit_t *hit) {
  assert(self != NULL);
  assert(seq != NULL);
  assert(seq_len > 0);
  assert(terminus_from + terminus_len <= seq_len);
  assert(read != NULL);
  assert(hit != NULL);
  /* (unpadded: the whole alignment is in the read, in_read not needed) */
  if (search(seq, seq_len, terminus_from, terminus_len, read, read_len, 0,
             (long)read_len, max_errors, hit, NULL))
    return true;
  if (!partial_left && !partial_right) {
    not_found(hit);
    return false;
  }
  size_t max_overhang = seq_len / 2;
  size_t pad_left = partial_left ? max_overhang : 0;
  size_t pad_right = partial_right ? max_overhang : 0;
  size_t padded_len = pad_left + read_len + pad_right;
  XCHECKSPACEINC(self->padded, &self->padded_cap, padded_len + 1, 1024);
  memset(self->padded, PAD_CHAR, pad_left);
  memcpy(self->padded + pad_left, read, read_len);
  memset(self->padded + pad_left + read_len, PAD_CHAR, pad_right);
  self->padded[padded_len] = '\0';
  long in_read;
  if (!search(seq, seq_len, terminus_from, terminus_len, self->padded,
              padded_len, (long)pad_left, (long)read_len, max_errors, hit,
              &in_read) ||
      hit->n_errors > floordiv(max_errors * in_read, (long)seq_len)) {
    not_found(hit);
    return false;
  }
  return true;
}

bool primer_search_run(primer_search_t *self, const primer_t *primer,
    const char *read, size_t read_len, long max_errors, bool partial_left,
    bool partial_right, primer_hit_t *hit) {
  assert(primer != NULL);
  return run(self, primer->seq, primer->seq_len, primer->terminus_from,
             primer->terminus_len, read, read_len, max_errors,
             partial_left, partial_right, hit);
}

bool primer_search_run_seq(primer_search_t *self, const char *seq,
    size_t seq_len, const char *read, size_t read_len, long max_errors,
    bool partial_left, bool partial_right, primer_hit_t *hit) {
  /* no terminus: the whole alignment is tail (see primer-search.h) */
  return run(self, seq, seq_len, seq_len, 0, read, read_len, max_errors,
             partial_left, partial_right, hit);
}

/* run a search and compare the result with the expected values
   (start == -1: not found) */
static void check(primer_search_t *s, const char *seq, const char *read,
    long max_errors, bool left, bool right, long start, long end,
    bool partial, long n_errors) {
  primer_hit_t hit;
  bool found = primer_search_run_seq(s, seq, strlen(seq), read,
      strlen(read), max_errors, left, right, &hit);
  test_assert(found == (start >= 0));
  test_eq(hit.start, start, long, "%ld");
  test_eq(hit.end, end, long, "%ld");
  if (found) {
    test_assert(hit.partial == partial);
    test_eq(hit.n_errors, n_errors, long, "%ld");
  }
}

/* the bounds of the terminus of a primer found in a read */
static void check_terminus(primer_search_t *s, const char *name,
    const char *read, long max_errors, bool left, bool right,
    long terminus_start, long terminus_end) {
  const primer_t *p = primers_builtin(name);
  test_assert(p != NULL);
  primer_hit_t hit;
  test_assert(primer_search_run(s, p, read, strlen(read), max_errors, left,
                                right, &hit));
  test_eq(hit.terminus_start, terminus_start, long, "%ld");
  test_eq(hit.terminus_end, terminus_end, long, "%ld");
}

void primer_search_test(void) {
  primer_search_t *s = primer_search_new();
  /* sequences without degenerate bases */
  check(s, "AGCTTA", "GGAGCTTAGG", 1, false, false, 2, 7, false, 0);
  check(s, "AGCTTA", "GGAGGTTAGG", 1, false, false, 2, 7, false, 1);
  check(s, "ACGTACGT", "", 2, false, false, -1, -1, false, 0);
  check(s, "ACGTACGT", "TACGT", 2, true, false, 0, 4, true, 0);
  check(s, "ACGTACGT", "GGACGTAC", 2, false, true, 2, 7, false, 2);
  check(s, "ACGTACGT", "GGGGGGACGTA", 2, false, true, 6, 10, true, 0);
  check(s, "ACGTACGT", "GGGGGGACGTA", 2, false, false, -1, -1, false, 0);
  /* degenerate bases of the primer (R = A or G) match either base */
  check(s, "ACRTACGT", "GGACATACGTGG", 0, false, false, 2, 9, false, 0);
  check(s, "ACRTACGT", "GGACGTACGTGG", 0, false, false, 2, 9, false, 0);
  check(s, "ACRTACGT", "GGACCTACGTGG", 0, false, false, -1, -1, false, 0);
  /* a degenerate base of the read (N) matches any base of the primer */
  check(s, "ACGTACGT", "GGACNTACGTGG", 0, false, false, 2, 9, false, 0);

  /*
     The terminus of a tailed primer: mbtuni-12 is
     ACGCGTGATC + AGCAAAAGCAGG (10 bases of tail, then the 12 bases of
     the conserved terminus, see primers.h)
  */
  /* the whole primer at the start of the read: the terminus begins at 10 */
  check_terminus(s, "mbtuni-12", "ACGCGTGATCAGCAAAAGCAGGTTTTTTTTTT", 2,
                 true, false, 10, 22);
  /* the primer after two bases of leader */
  check_terminus(s, "mbtuni-12", "TTACGCGTGATCAGCAAAAGCAGGTTTTTTTT", 2,
                 true, false, 12, 24);
  /* a deletion in the tail moves the terminus one base to the left */
  check_terminus(s, "mbtuni-12", "ACGCGTGAT" "AGCAAAAGCAGGTTTTTTTTTT", 2,
                 true, false, 9, 21);
  /* an insertion in the tail moves it one base to the right */
  check_terminus(s, "mbtuni-12", "ACGCGTGATCCAGCAAAAGCAGGTTTTTTTTT", 2,
                 true, false, 11, 23);
  /* the read begins inside the terminus (the primer overhangs the read
     on the left, here by half of its length, the maximum): the terminus
     begins one base before the read */
  check_terminus(s, "mbtuni-12", "GCAAAAGCAGGTTTTTTTTTTTTTTTTTTTT", 2,
                 true, false, -1, 11);
  /* a longer overhang (only 9 of the 22 bases of the primer are in the
     read) is not accepted */
  {
    const primer_t *p = primers_builtin("mbtuni-12");
    primer_hit_t hit;
    test_assert(!primer_search_run(s, p, "AAAAGCAGGTTTTTTTTTTTTTTTTTTTT",
                                   29, 2, true, false, &hit));
  }
  /* a primer without a tail: the terminus is the whole alignment */
  check_terminus(s, "uni12", "TTAGCAAAAGCAGGTTTTTTTTTTTTTT", 2, true,
                 false, 2, 14);
  /* the reverse complement of a primer, as it is found at the other end
     of an amplicon: its terminus is its prefix */
  {
    const primer_t *p = primers_builtin("mbtuni-13");
    char rc_seq[32];
    primer_t rc;
    primers_reverse(p, rc_seq, &rc);
    test_eq_size(rc.terminus_from, 0);
    test_eq_mem(rc.seq, rc.seq_len, "CCTTGTTTCTACTGATCACGCGT");
    primer_hit_t hit;
    /* 10 bases of segment, the terminus, then the reverse tail */
    const char *read = "TTTTTTTTTT" "CCTTGTTTCTACT" "GATCACGCGT";
    test_assert(primer_search_run(s, &rc, read, strlen(read), 2, false,
                                  true, &hit));
    test_eq(hit.start, 10, long, "%ld");
    test_eq(hit.terminus_start, 10, long, "%ld");
    test_eq(hit.terminus_end, 23, long, "%ld");
  }
  primer_search_delete(s);
}
