/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "fasta.h"
#include "junction.h"
#include "test.h"

bool junction_expected(long len, long distance) {
  return distance - len > JUNCTION_MIN_DELETION + distance / 10;
}

/* the flank to which a probe belongs */
typedef enum { FLANK_NONE, FLANK_LEFT, FLANK_RIGHT } flank_t;

/* the search of the junction of a template */
typedef struct {
  reference_set_t *refs;
  const char *seq;
  long len;
  size_t ref_index;
  long probe_len;
  /* the diagonals of the flanks, as given by the ends of the template:
     the position on the reference of the base 0 of the template, if the
     template continued along that flank */
  long left_diagonal, right_diagonal;
  /* a probe belongs to a flank if its diagonal is within this distance
     of that of the flank (the indels of the read move it) */
  long tolerance;
} search_t;

/* the flank of the probe at <offset> of the template; <diagonal>: that
   of the probe, if it belongs to a flank */
static flank_t probe_flank(const search_t *s, long offset, long *diagonal) {
  reference_hit_t hit;
  reference_set_locate_in(s->refs, s->seq + offset,
      (size_t)(s->len - offset), s->ref_index, false, &hit);
  if (!hit.mapped)
    return FLANK_NONE;
  long d = hit.ref_start - offset;
  long from_left = d - s->left_diagonal, from_right = d - s->right_diagonal;
  if (from_left < 0) from_left = -from_left;
  if (from_right < 0) from_right = -from_right;
  if (from_left > s->tolerance && from_right > s->tolerance)
    return FLANK_NONE;
  *diagonal = d;
  return from_left <= from_right ? FLANK_LEFT : FLANK_RIGHT;
}

/* the greatest share of the bases of a template which may differ from
   the reference in its two flanks, for a junction found by alignment */
#define JUNCTION_MAX_DIFFERENCE 0.3

/*
   Align bases of the template to the reference from one end of them,
   within a band of <band> bases: the <n> bases from <from> on, forwards
   from the first one at <anchor> of the reference, or (<backwards>) the
   <n> bases before <from> + <n>, from the last one at <anchor> towards
   the start of both. For each number <i> of bases, from 0 to <n>:
   cost[i], the fewest differences (edit distance) of those bases to
   bases of the reference from the anchor on, and used[i], how many
   bases of the reference that takes.
*/
static void flank_costs(const search_t *s, long from, long n, long anchor,
    bool backwards, long band, long *cost, long *used) {
  const fasta_record_t *ref =
      reference_set_fasta(s->refs)->records + s->ref_index;
  const long w = 2 * band + 1, inf = n + w;
  long *prev, *row;
  XMALLOC(prev, sizeof(long) * 2 * (size_t)w);
  row = prev + w;
  /* (d: bases of the reference less bases of the template, plus band) */
  for (long d = 0; d < w; d++)
    prev[d] = d >= band ? d - band : inf;
  cost[0] = 0;
  used[0] = 0;
  for (long i = 1; i <= n; i++) {
    char base = backwards ? s->seq[from + n - i] : s->seq[from + i - 1];
    cost[i] = inf;
    for (long d = 0; d < w; d++) {
      long j = i + d - band;   /* bases of the reference */
      long best = inf;
      if (j >= 0) {
        if (j >= 1) {
          long p = backwards ? anchor - (j - 1) : anchor + (j - 1);
          bool differs = p < 0 || (size_t)p >= ref->seq_len ||
                         ref->seq[p] != base;
          best = prev[d] + differs;
          if (d > 0 && row[d - 1] + 1 < best)
            best = row[d - 1] + 1;
        }
        if (d + 1 < w && prev[d + 1] + 1 < best)
          best = prev[d + 1] + 1;
      }
      row[d] = best;
      if (best < cost[i]) {
        cost[i] = best;
        used[i] = j;
      }
    }
    memcpy(prev, row, sizeof(long) * (size_t)w);
  }
  free(prev);
}

/* locate the junction among the bases <from> to <to> - 1 of the
   template by aligning them from both sides: the bases before the
   junction along the diagonal <left>, the others along <right> (see
   search_t). The junction is where the differences of the two parts
   are the fewest; of equally good positions, the first (the bases which
   fit both flanks are given to the right one). The band of the
   alignments grows with the number of bases (the indels of the read) */
static void locate_by_alignment(const search_t *s, long from, long to,
    long left, long right, junction_t *out) {
  long n = to - from, band = 8 + n / 8;
  if (n <= 0)
    return;
  long *buf;
  XMALLOC(buf, sizeof(long) * 4 * (size_t)(n + 1));
  long *lcost = buf, *lused = buf + (n + 1), *rcost = buf + 2 * (n + 1),
       *rused = buf + 3 * (n + 1);
  flank_costs(s, from, n, left + from, false, band, lcost, lused);
  flank_costs(s, from, n, right + to - 1, true, band, rcost, rused);
  long split = 0;
  for (long i = 1; i <= n; i++)
    if (lcost[i] + rcost[n - i] < lcost[split] + rcost[n - split])
      split = i;
  long left_to = left + from + lused[split],
       right_from = right + to - rused[n - split];
  if ((double)(lcost[split] + rcost[n - split]) <=
          JUNCTION_MAX_DIFFERENCE * (double)n &&
      right_from - left_to >= JUNCTION_MIN_DELETION) {
    out->found = true;
    out->left_to = left_to;
    out->right_from = right_from;
  }
  free(buf);
}

long junction_flank_differences(reference_set_t *refs, const char *seq,
    size_t len, size_t ref_index, long anchor, bool backwards) {
  assert(refs != NULL);
  assert(seq != NULL || len == 0);
  if (len == 0)
    return 0;
  search_t s = {.refs = refs, .seq = seq, .len = (long)len,
                .ref_index = ref_index};
  long *buf;
  XMALLOC(buf, sizeof(long) * 2 * (len + 1));
  flank_costs(&s, 0, (long)len, anchor, backwards, 8 + (long)len / 8, buf,
              buf + len + 1);
  long differences = buf[len];
  free(buf);
  return differences;
}

/*
   Locate the junction (see junction.h). Probes tell on which flank each
   part of the template lies, by bisection, as long as they can: the
   junction is then between the last probe of the left flank and the
   first of the right one, and is located among those bases by alignment
   (locate_by_alignment()). A probe which does not tell its flank (it
   contains the junction, or the sample differs too much from the
   reference there) ends the bisection: the alignment is then over the
   bases between the probes known so far, at first the whole template.
*/
void junction_locate(reference_set_t *refs, const char *seq, size_t len,
    size_t ref_index, long first, long last, junction_t *out) {
  assert(refs != NULL);
  assert(seq != NULL || len == 0);
  assert(out != NULL);
  out->found = false;
  long n = (long)len, deleted = last - first + 1 - n;
  if (n == 0 || !junction_expected(n, last - first + 1))
    return;
  search_t s = {.refs = refs, .seq = seq, .len = n, .ref_index = ref_index,
                .probe_len = (long)reference_set_probe_len(refs),
                .left_diagonal = first, .right_diagonal = last - (n - 1)};
  s.tolerance = 20 + n / 10;
  if (s.tolerance > (deleted - 1) / 2)
    s.tolerance = (deleted - 1) / 2;
  /* the last probe known to be of the left flank and the first known to
     be of the right flank, with their diagonals: at first the two ends
     of the template */
  long lo = 0, hi = n > s.probe_len ? n - s.probe_len : 0;
  long left = s.left_diagonal, right = s.right_diagonal;
  bool told = true;
  while (told && hi - lo > 2 * s.probe_len + 1) {
    /* the probe in the middle, or, if it does not tell, the one after
       it or before it: all three are between the two known probes,
       which thus get closer at each step */
    long mid = lo + (hi - lo) / 2;
    const long tries[3] = {mid, mid + s.probe_len, mid - s.probe_len};
    flank_t flank = FLANK_NONE;
    for (size_t i = 0; i < 3 && flank == FLANK_NONE; i++) {
      long diagonal;
      flank = probe_flank(&s, tries[i], &diagonal);
      if (flank == FLANK_LEFT) {
        lo = tries[i];
        left = diagonal;
      } else if (flank == FLANK_RIGHT) {
        hi = tries[i];
        right = diagonal;
      }
    }
    told = flank != FLANK_NONE;
  }
  long to = hi + s.probe_len < n ? hi + s.probe_len : n;
  locate_by_alignment(&s, lo, to, left, right, out);
}

/* --- unit tests --- */

/* a deterministic pseudo-random sequence of <len> bases into <seq> */
static void random_seq(char *seq, size_t len, unsigned long *state) {
  for (size_t i = 0; i < len; i++) {
    *state = *state * 6364136223846793005UL + 1442695040888963407UL;
    seq[i] = "ACGT"[(*state >> 33) & 3];
  }
  seq[len] = '\0';
}

/* the template made of the first <left> bases of <ref> and of its bases
   from <right> on, into <read>; its length */
static size_t deleted_copy(const char *ref, size_t ref_len, size_t left,
    size_t right, char *read) {
  memcpy(read, ref, left);
  memcpy(read + left, ref + right, ref_len - right);
  return left + ref_len - right;
}

void junction_test(void) {
  enum { REF_LEN = 1000, PROBE = 32 };
  unsigned long state = 11;
  char ref0[REF_LEN + 1], ref1[REF_LEN + 1], read[REF_LEN + 1];
  random_seq(ref0, REF_LEN, &state);
  random_seq(ref1, REF_LEN, &state);
  /* the last base before each deletion tested below differs from the
     last base deleted: the junction is then at one position only */
  ref1[149] = 'A'; ref1[599] = 'C';
  ref1[19] = 'A'; ref1[499] = 'C';
  ref1[399] = 'A'; ref1[979] = 'C';
  ref1[299] = 'A'; ref1[699] = 'C';
  fasta_record_t recs[2] = {{"r0", "", ref0, REF_LEN},
                            {"r1", "", ref1, REF_LEN}};
  fasta_t fa = {recs, 2};
  reference_index_t *ix = reference_index_new(&fa);
  reference_set_t *rs = reference_set_new(ix, PROBE, 0.2);
  junction_t j;

  test_assert(junction_expected(550, 1000));
  test_assert(!junction_expected(950, 1000));   /* indels of a long read */
  test_assert(!junction_expected(1000, 1000));

  /* the bases 150 to 599 of the reference are missing */
  size_t len = deleted_copy(ref1, REF_LEN, 150, 600, read);
  junction_locate(rs, read, len, 1, 0, REF_LEN - 1, &j);
  test_assert(j.found);
  test_eq_long(j.left_to, 150);
  test_eq_long(j.right_from, 600);

  /* the same, with errors in the read on both sides of the junction */
  read[100] = read[100] == 'A' ? 'C' : 'A';
  read[200] = read[200] == 'A' ? 'C' : 'A';
  junction_locate(rs, read, len, 1, 0, REF_LEN - 1, &j);
  test_assert(j.found);
  test_eq_long(j.left_to, 150);
  test_eq_long(j.right_from, 600);

  /* a flank shorter than a probe, at either end */
  len = deleted_copy(ref1, REF_LEN, 20, 500, read);
  junction_locate(rs, read, len, 1, 0, REF_LEN - 1, &j);
  test_assert(j.found);
  test_eq_long(j.left_to, 20);
  test_eq_long(j.right_from, 500);
  len = deleted_copy(ref1, REF_LEN, 400, 980, read);
  junction_locate(rs, read, len, 1, 0, REF_LEN - 1, &j);
  test_assert(j.found);
  test_eq_long(j.left_to, 400);
  test_eq_long(j.right_from, 980);

  /* a template which is a part of the reference: 100 to 299 and 700 to
     899 */
  memcpy(read, ref1 + 100, 200);
  memcpy(read + 200, ref1 + 700, 200);
  junction_locate(rs, read, 400, 1, 100, 899, &j);
  test_assert(j.found);
  test_eq_long(j.left_to, 300);
  test_eq_long(j.right_from, 700);

  /* deletions of 400 bases at many positions (templates of many
     lengths): the two flanks found are the template (the junction can be
     a few bases before where it was made, if the bases before it are the
     same on its two sides) */
  char joined[REF_LEN + 1];
  size_t n_wrong = 0;
  for (size_t left = 33; left < 300; left++) {
    len = deleted_copy(ref1, REF_LEN, left, left + 400, read);
    junction_locate(rs, read, len, 1, 0, REF_LEN - 1, &j);
    if (!j.found || j.left_to < 0 || j.right_from - j.left_to != 400 ||
        deleted_copy(ref1, REF_LEN, (size_t)j.left_to,
                     (size_t)j.right_from, joined) != len ||
        memcmp(joined, read, len) != 0)
      n_wrong++;
  }
  test_eq_size(n_wrong, 0);

  /* no deletion */
  junction_locate(rs, ref1, REF_LEN, 1, 0, REF_LEN - 1, &j);
  test_assert(!j.found);

  /* what lies between the flanks is not in the reference: no probe
     tells where the junction is */
  len = deleted_copy(ref1, REF_LEN, 150, 600, read);
  random_seq(read + 100, 350, &state);
  read[450] = ref1[900];   /* (random_seq() ends the sequence) */
  junction_locate(rs, read, len, 1, 0, REF_LEN - 1, &j);
  test_assert(!j.found);

  reference_set_delete(rs);
  reference_index_delete(ix);
}
