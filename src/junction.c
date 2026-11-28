/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "fasta.h"
#include "junction.h"
#include "test.h"

bool junction_expected(long len, long distance) {
  long allowance = distance / 10;
  if (allowance < JUNCTION_MIN_DELETION)
    allowance = JUNCTION_MIN_DELETION;
  return distance - len > allowance;
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

/* locate the probe at <offset> of the template where it lies if it
   continues the flank of diagonal <flank>; false if it is not there */
static bool probe_on_flank(const search_t *s, long offset, long flank,
    reference_hit_t *hit) {
  size_t len = (size_t)(s->len - offset);
  if (len > (size_t)s->probe_len)
    len = (size_t)s->probe_len;
  /* (the first base of the probe within the tolerance of its expected
     position, and room for the rest of it, with its insertions) */
  long expected = flank + offset;
  reference_set_locate_within(s->refs, s->seq + offset, len, s->ref_index,
      false, expected - s->tolerance,
      expected + s->tolerance + 2 * (long)len, hit);
  if (!hit->mapped)
    return false;
  long distance = hit->ref_start - expected;
  return distance >= -s->tolerance && distance <= s->tolerance;
}

/* the flank of the probe at <offset> of the template; <diagonal>: that
   of the probe, if it belongs to a flank. The probe is aligned only
   where the two flanks would have it: of two matches, the one with
   fewer errors, then that of the left flank */
static flank_t probe_flank(const search_t *s, long offset, long *diagonal) {
  reference_hit_t left, right;
  bool on_left = probe_on_flank(s, offset, s->left_diagonal, &left);
  bool on_right = probe_on_flank(s, offset, s->right_diagonal, &right);
  if (!on_left && !on_right)
    return FLANK_NONE;
  if (on_left && (!on_right || left.n_errors <= right.n_errors)) {
    *diagonal = left.ref_start - offset;
    return FLANK_LEFT;
  }
  *diagonal = right.ref_start - offset;
  return FLANK_RIGHT;
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
  const long w = 2 * band + 1, max_ref = n + band;
  const int inf = (int)(n + w);
  /* the bases of the reference, in the order in which they are aligned
     (0 where there is none: it differs from any base), and two rows of
     the matrix, with a cell before and after each (inf: no path) */
  char *bases;
  int *rows;
  XMALLOC(bases, (size_t)max_ref);
  XMALLOC(rows, sizeof(int) * 2 * (size_t)(w + 2));
  for (long k = 0; k < max_ref; k++) {
    long p = backwards ? anchor - k : anchor + k;
    bases[k] = p < 0 || (size_t)p >= ref->seq_len ? 0 : ref->seq[p];
  }
  int *prev = rows + 1, *row = rows + (w + 2) + 1;
  prev[-1] = prev[w] = row[-1] = row[w] = inf;
  /* (d: bases of the reference less bases of the template, plus band) */
  for (long d = 0; d < w; d++)
    prev[d] = d >= band ? (int)(d - band) : inf;
  cost[0] = 0;
  used[0] = 0;
  for (long i = 1; i <= n; i++) {
    char base = backwards ? s->seq[from + n - i] : s->seq[from + i - 1];
    /* (j = i + d - band bases of the reference: none below <first>) */
    long first = band > i ? band - i : 0;
    int lowest = inf;
    long lowest_d = first;
    for (long d = 0; d < first; d++)
      row[d] = inf;
    if (i <= band) {
      /* no base of the reference */
      row[first] = prev[first + 1] + 1;
      lowest = row[first];
      first++;
    }
    /* (the base of the reference aligned in the cell d) */
    const char *aligned = bases + (i - band - 1);
    for (long d = first; d < w; d++) {
      int best = prev[d] + (aligned[d] != base);
      if (row[d - 1] + 1 < best)
        best = row[d - 1] + 1;
      if (prev[d + 1] + 1 < best)
        best = prev[d + 1] + 1;
      row[d] = best;
      if (best < lowest) {
        lowest = best;
        lowest_d = d;
      }
    }
    cost[i] = lowest;
    used[i] = i + lowest_d - band;
    int *swap = prev;
    prev = row;
    row = swap;
  }
  free(bases);
  free(rows);
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
    out->offset = from + split;
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
  test_assert(junction_expected(850, 1000));    /* a short deletion */
  test_assert(!junction_expected(900, 1000));   /* the smallest deletion */
  test_assert(!junction_expected(950, 1000));
  test_assert(!junction_expected(1000, 1000));
  test_assert(!junction_expected(2100, 2300));  /* indels of a long read */
  test_assert(junction_expected(2000, 2300));

  /* the bases 150 to 599 of the reference are missing */
  size_t len = deleted_copy(ref1, REF_LEN, 150, 600, read);
  junction_locate(rs, read, len, 1, 0, REF_LEN - 1, &j);
  test_assert(j.found);
  test_eq_long(j.left_to, 150);
  test_eq_long(j.right_from, 600);
  test_eq_long(j.offset, 150);

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
  test_eq_long(j.offset, 200);

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
