/* (c) 2026 Giorgio Gonnella */

/*
   Unit tests of the scan of a sequence for a pattern (pattern-scan.h):
   examples with a primer terminus, and random patterns and sequences
   compared with the plain dynamic programming.
*/

#include <assert.h>
#include <stdbool.h>
#include <string.h>
#include "pattern-scan.h"
#include "pattern-scan-test.h"
#include "iupac.h"
#include "test.h"

/* whether the characters <a> and <b> are equal (see iupac_equalities) */
static bool chars_equal(char a, char b) {
  if (a == b)
    return true;
  for (int i = 0; i < n_iupac_equalities; i++) {
    const EdlibEqualityPair *e = iupac_equalities + i;
    if ((e->first == a && e->second == b) ||
        (e->first == b && e->second == a))
      return true;
  }
  return false;
}

/* the edit distance of <pattern> (length <m>) to the bases of <seq>
   which end at each position (into <out>), by the plain dynamic
   programming with chars_equal() */
static void distances(const char *pattern, size_t m, const char *seq,
    size_t n, long *out) {
  long col[2 * PATTERN_SCAN_WORD + 1];
  for (size_t i = 0; i <= m; i++)
    col[i] = (long)i;
  for (size_t j = 0; j < n; j++) {
    long diagonal = col[0];
    for (size_t i = 1; i <= m; i++) {
      long best = diagonal + (chars_equal(pattern[i - 1], seq[j]) ? 0 : 1);
      if (col[i] + 1 < best) best = col[i] + 1;
      if (col[i - 1] + 1 < best) best = col[i - 1] + 1;
      diagonal = col[i];
      col[i] = best;
    }
    out[j] = col[m];
  }
}

/* the next number of a random sequence */
static unsigned long next_random(unsigned long *state) {
  *state = *state * 6364136223846793005UL + 1442695040888963407UL;
  return *state >> 33;
}

/* fill <seq> with <len> random characters of <alphabet> */
static void random_chars(char *seq, size_t len, const char *alphabet,
    unsigned long *state) {
  size_t n = strlen(alphabet);
  for (size_t i = 0; i < len; i++)
    seq[i] = alphabet[next_random(state) % n];
  seq[len] = '\0';
}

/* whether a run of end positions and a match are as given */
static bool end_is(const pattern_scan_end_t *e, long first, long best,
    long last, long n_errors) {
  return e->first == first && e->best == best && e->last == last &&
         e->n_errors == n_errors;
}

static bool match_is(const pattern_match_t *m, long start, long end,
    long n_errors) {
  return m->start == start && m->end == end && m->n_errors == n_errors;
}

/* scan <seq> for <pattern> and compare the runs with the distances
   computed by distances() */
static void check_scan(const char *pattern, const char *seq, long k) {
  enum { MAX = 512 };
  size_t m = strlen(pattern), n = strlen(seq);
  assert(n <= MAX);
  long expected[MAX];
  distances(pattern, m, seq, n, expected);
  pattern_scan_t *p = pattern_scan_new(pattern, m);
  pattern_scan_end_t ends[MAX];
  size_t n_ends = pattern_scan_run(p, seq, n, k, ends, MAX);
  size_t e = 0;
  bool same = true;
  for (size_t j = 0; j < n;) {
    if (expected[j] > k) {
      j++;
      continue;
    }
    size_t last = j, best = j;
    while (last + 1 < n && expected[last + 1] <= k)
      last++;
    for (size_t i = j; i <= last; i++)
      if (expected[i] <= expected[best])
        best = i;
    same = same && e < n_ends && ends[e].first == (long)j &&
           ends[e].last == (long)last && ends[e].best == (long)best &&
           ends[e].n_errors == expected[best];
    e++;
    j = last + 1;
  }
  test_assert(same && e == n_ends);
  pattern_scan_delete(p);
}

/* scan <seq> for <back> (backwards) and <forth> (forwards) in one pass
   and compare with the two scans done one after the other */
static void check_both(const char *back, const char *forth,
    const char *seq, long k) {
  size_t n = strlen(seq);
  pattern_scan_t *pb = pattern_scan_new(back, strlen(back));
  pattern_scan_t *pf = pattern_scan_new(forth, strlen(forth));
  pattern_scan_best_t b, f, eb = {0}, ef = {0};
  eb.found = pattern_scan_best(pb, seq, n, k, true, &eb.n_errors,
                               &eb.first, &eb.last);
  ef.found = pattern_scan_best(pf, seq, n, k, false, &ef.n_errors,
                               &ef.first, &ef.last);
  pattern_scan_best_both(pb, pf, seq, n, k, &b, &f);
  test_assert(b.found == eb.found && f.found == ef.found &&
      (!b.found || (b.n_errors == eb.n_errors && b.first == eb.first &&
                    b.last == eb.last)) &&
      (!f.found || (f.n_errors == ef.n_errors && f.first == ef.first &&
                    f.last == ef.last)));
  pattern_scan_delete(pb);
  pattern_scan_delete(pf);
}

/* two patterns in one pass, against one pass for each: patterns taken
   from the sequence (found) and random ones (found or not), of
   different lengths, also longer than a word */
static void both_test(void) {
  unsigned long state = 11;
  char back[2 * PATTERN_SCAN_WORD + 1], forth[2 * PATTERN_SCAN_WORD + 1];
  char text[301];
  const size_t lens[] = {1, 8, 32, 64, 65};
  const size_t n_lens = sizeof(lens) / sizeof(lens[0]);
  for (int rep = 0; rep < 200; rep++) {
    size_t n = rep == 0 ? 0 : 150 + (size_t)(rep % 150);
    random_chars(text, n, rep % 3 == 0 ? "AC" : "ACGT", &state);
    size_t mb = lens[rep % n_lens], mf = lens[(rep / n_lens) % n_lens];
    random_chars(back, mb, "ACGTN", &state);
    random_chars(forth, mf, "ACGTR", &state);
    /* (the backwards pattern is given reversed) */
    if (rep % 2 == 0 && n > 0)
      for (size_t i = 0; i < mb; i++)
        back[i] = text[70 + mb - 1 - i];
    if (rep % 4 < 2 && n > 0)
      memcpy(forth, text + 20, mf);
    for (long k = 0; k <= 6; k += 3)
      check_both(back, forth, text, k);
  }
}

void pattern_scan_test(void) {
  both_test();
  const char *uni12 = "AGCRAAAGCAGG";
  pattern_scan_t *p = pattern_scan_new(uni12, strlen(uni12));
  pattern_scan_end_t ends[4];
  /* exact, with either base of the degenerate position */
  const char *seq = "TTTTAGCAAAAGCAGGTTTTTTTTAGCGAAAGCAGGTTTT";
  test_eq_size(pattern_scan_run(p, seq, strlen(seq), 0, ends, 4), 2);
  test_assert(end_is(ends, 15, 15, 15, 0));
  test_assert(end_is(ends + 1, 35, 35, 35, 0));
  /* with one error allowed, the positions next to the end are
     acceptable too: one run for each match */
  test_eq_size(pattern_scan_run(p, seq, strlen(seq), 1, ends, 4), 2);
  test_assert(end_is(ends, 14, 15, 16, 0));
  /* no room for the second run: not reported */
  test_eq_size(pattern_scan_run(p, seq, strlen(seq), 1, ends, 1), 1);
  test_eq_long(ends[0].best, 15);
  /* a mismatch, a deletion, an N in the sequence; nothing in a sequence
     without the pattern, or in an empty one */
  seq = "CCAGCAATAGCAGGCCCCAGCAAAGCAGGCCCCAGCANAAGCAGG";
  test_eq_size(pattern_scan_run(p, seq, strlen(seq), 1, ends, 4), 3);
  test_assert(ends[0].n_errors == 1 && ends[1].n_errors == 1 &&
              ends[2].n_errors == 0);
  test_eq_size(pattern_scan_run(p, "CCCCCCCCCCCCCCCCCCCC", 20, 2, ends, 4),
               0);
  test_eq_size(pattern_scan_run(p, "", 0, 2, ends, 4), 0);
  /* the matches themselves: where they begin as well */
  pattern_match_t matches[PATTERN_SCAN_MAX_MATCHES];
  test_eq_size(pattern_scan_matches(p, seq, strlen(seq), 1, matches), 3);
  test_assert(match_is(matches, 2, 13, 1));
  test_assert(match_is(matches + 1, 18, 28, 1));
  test_assert(match_is(matches + 2, 33, 44, 0));
  pattern_scan_delete(p);

  /* against the plain dynamic programming: random sequences with the
     pattern planted with errors, for patterns of one base, of a word
     and longer than a word */
  unsigned long state = 7;
  char pattern[2 * PATTERN_SCAN_WORD + 1], text[401];
  const size_t lens[] = {1, 12, 23, 63, 64, 65, 100};
  for (size_t l = 0; l < sizeof(lens) / sizeof(lens[0]); l++)
    for (int rep = 0; rep < 20; rep++) {
      size_t m = lens[l];
      random_chars(pattern, m, "ACGTRN", &state);
      size_t n = 0;
      while (n + m + 20 < sizeof(text)) {
        for (size_t i = 0; i < 20 + m; i++) {
          unsigned long r = next_random(&state);
          /* random bases, then the pattern with some bases changed */
          text[n++] = i < 20 || r % 10 == 0 || pattern[i - 20] == 'R' ?
                      "ACGT"[r % 4] : pattern[i - 20];
        }
      }
      text[n] = '\0';
      for (long k = 0; k <= 3 && k < (long)m; k++)
        check_scan(pattern, text, k);
    }
}
