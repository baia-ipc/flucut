/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "pattern-scan.h"
#include "primer-search.h"
#include "test.h"

struct pattern_scan_t {
  char *seq;
  size_t len, cap;
  /* for each character, the positions of the pattern it is equal to:
     as the bits of a word (a pattern of at most PATTERN_SCAN_WORD
     bases), or one byte for each position (a longer one) */
  uint64_t peq[256];
  unsigned char *equal;
  /* a column of the matrix of the edit distances (len + 1) */
  long *col;
};

/* for each character, the characters equal to it (see iupac_equalities):
   computed once, shared by all the patterns */
#define MAX_EQUAL 24
static unsigned char equal_chars[256][MAX_EQUAL];
static unsigned char n_equal_chars[256];
static pthread_once_t equal_chars_once = PTHREAD_ONCE_INIT;

static void init_equal_chars(void) {
  for (int c = 0; c < 256; c++) {
    equal_chars[c][0] = (unsigned char)c;
    n_equal_chars[c] = 1;
  }
  for (int i = 0; i < n_iupac_equalities; i++) {
    unsigned char pair[2] = {(unsigned char)iupac_equalities[i].first,
                             (unsigned char)iupac_equalities[i].second};
    for (int side = 0; side < 2; side++) {
      unsigned char a = pair[side], b = pair[1 - side];
      bool known = false;
      for (int j = 0; j < n_equal_chars[a]; j++)
        known = known || equal_chars[a][j] == b;
      if (!known) {
        assert(n_equal_chars[a] < MAX_EQUAL);
        equal_chars[a][n_equal_chars[a]++] = b;
      }
    }
  }
}

void pattern_scan_set(pattern_scan_t *self, const char *pattern,
    size_t len) {
  assert(self != NULL);
  assert(pattern != NULL);
  assert(len > 0);
  if (len > self->cap) {
    self->cap = len;
    XREALLOC(self->seq, len + 1);
    XREALLOC(self->col, sizeof(long) * (len + 1));
    if (len > PATTERN_SCAN_WORD)
      XREALLOC(self->equal, 256 * len);
  }
  memcpy(self->seq, pattern, len);
  self->seq[len] = '\0';
  self->len = len;
  bool word = len <= PATTERN_SCAN_WORD;
  if (word)
    memset(self->peq, 0, sizeof(self->peq));
  else
    memset(self->equal, 0, 256 * len);
  for (size_t i = 0; i < len; i++) {
    unsigned char a = (unsigned char)pattern[i];
    for (int j = 0; j < n_equal_chars[a]; j++) {
      unsigned char c = equal_chars[a][j];
      if (word)
        self->peq[c] |= (uint64_t)1 << i;
      else
        self->equal[(size_t)c * len + i] = 1;
    }
  }
}

pattern_scan_t* pattern_scan_new(const char *pattern, size_t len) {
  pthread_once(&equal_chars_once, init_equal_chars);
  pattern_scan_t *self;
  XMALLOC(self, sizeof(pattern_scan_t));
  self->seq = NULL;
  self->equal = NULL;
  self->col = NULL;
  self->cap = 0;
  pattern_scan_set(self, pattern, len);
  return self;
}

void pattern_scan_delete(pattern_scan_t *self) {
  if (self == NULL)
    return;
  free(self->equal);
  free(self->col);
  free(self->seq);
  free(self);
}

/* whether the base <i> of the pattern is equal to the character <c> */
static inline bool equal_at(const pattern_scan_t *self, size_t i,
    unsigned char c) {
  return self->len <= PATTERN_SCAN_WORD ? (self->peq[c] >> i) & 1 :
         self->equal[(size_t)c * self->len + i];
}

/*
   The edit distance of the pattern to the bases which end at each
   position, one position at a time.

   With one word (Myers 1999, as formulated by Hyyro 2001 for the search
   of a pattern: the alignment may begin anywhere in the sequence):
   <pv> and <mv> are the positions of the pattern at which the distance
   grows, or decreases, by one from the previous position of the
   pattern, for the bases read so far. A longer pattern: the column of
   the matrix of the distances (self->col).
*/
typedef struct {
  uint64_t pv, mv;
  long score;
} scan_state_t;

static void scan_start(pattern_scan_t *self, scan_state_t *st) {
  st->pv = ~(uint64_t)0;
  st->mv = 0;
  st->score = (long)self->len;
  if (self->len > PATTERN_SCAN_WORD)
    for (size_t i = 0; i <= self->len; i++)
      self->col[i] = (long)i;
}

/* a step of the scan with one word: <eq>, the positions of the pattern
   equal to the character read; into <up> and <down>, the positions of
   the pattern at which the distance has grown, or decreased, by one
   from the previous base of the sequence */
static inline void word_step(uint64_t eq, uint64_t *pv, uint64_t *mv,
    uint64_t *up, uint64_t *down) {
  uint64_t xv = eq | *mv;
  uint64_t xh = (((eq & *pv) + *pv) ^ *pv) | eq;
  uint64_t ph = *mv | ~(xh | *pv);
  uint64_t mh = *pv & xh;
  *up = ph;
  *down = mh;
  ph <<= 1;
  mh <<= 1;
  *pv = mh | ~(xv | ph);
  *mv = ph & xv;
}

/* read the character <c>: the distance of the pattern to the bases
   which end with it */
static inline long scan_step(pattern_scan_t *self, scan_state_t *st,
    unsigned char c) {
  const size_t m = self->len;
  if (m <= PATTERN_SCAN_WORD) {
    uint64_t up, down;
    word_step(self->peq[c], &st->pv, &st->mv, &up, &down);
    /* (without a branch: the distance changes at about every other
       base) */
    st->score += (long)((up >> (m - 1)) & 1) -
                 (long)((down >> (m - 1)) & 1);
    return st->score;
  }
  const unsigned char *equal = self->equal + (size_t)c * m;
  long *col = self->col;
  long diagonal = col[0];   /* (col[0] stays 0: a free start) */
  for (size_t i = 1; i <= m; i++) {
    long best = diagonal + (equal[i - 1] ? 0 : 1);
    if (col[i] + 1 < best)
      best = col[i] + 1;
    if (col[i - 1] + 1 < best)
      best = col[i - 1] + 1;
    diagonal = col[i];
    col[i] = best;
  }
  return col[m];
}

size_t pattern_scan_run(pattern_scan_t *self, const char *seq, size_t len,
    long max_errors, pattern_scan_end_t *ends, size_t max_ends) {
  assert(self != NULL);
  assert(seq != NULL || len == 0);
  assert(max_errors >= 0 && max_errors < (long)self->len);
  scan_state_t st;
  scan_start(self, &st);
  size_t n = 0;
  bool open = false;
  for (size_t j = 0; j < len; j++) {
    long score = scan_step(self, &st, (unsigned char)seq[j]);
    if (score > max_errors) {
      open = false;
      continue;
    }
    if (!open) {
      open = true;
      /* (no room: the run is not reported) */
      if (n++ >= max_ends)
        continue;
      ends[n - 1].first = (long)j;
      ends[n - 1].n_errors = score + 1;
    }
    if (n > max_ends)
      continue;
    pattern_scan_end_t *e = ends + n - 1;
    e->last = (long)j;
    if (score <= e->n_errors) {
      e->best = (long)j;
      e->n_errors = score;
    }
  }
  return n < max_ends ? n : max_ends;
}

bool pattern_scan_best(pattern_scan_t *self, const char *seq, size_t len,
    long max_errors, bool backwards, long *n_errors, long *first,
    long *last) {
  assert(self != NULL);
  assert(seq != NULL || len == 0);
  assert(max_errors >= 0);
  scan_state_t st;
  scan_start(self, &st);
  long best = max_errors + 1;
  for (size_t t = 0; t < len; t++) {
    size_t j = backwards ? len - 1 - t : t;
    long score = scan_step(self, &st, (unsigned char)seq[j]);
    if (score > best)
      continue;
    if (score < best) {
      best = score;
      *first = *last = (long)j;
    } else if (score <= max_errors) {
      if (backwards)
        *first = (long)j;
      else
        *last = (long)j;
    }
  }
  if (best > max_errors)
    return false;
  *n_errors = best;
  return true;
}

/* note the distance <score> at the position <j> (see
   pattern_scan_best()); <lowest>: the lowest so far */
static inline void note_score(long score, long j, long max_errors,
    bool backwards, long *lowest, pattern_scan_best_t *best) {
  if (score > *lowest)
    return;
  if (score < *lowest) {
    *lowest = score;
    best->first = best->last = j;
  } else if (score <= max_errors) {
    if (backwards)
      best->first = j;
    else
      best->last = j;
  }
}

void pattern_scan_best_both(pattern_scan_t *backwards,
    pattern_scan_t *forwards, const char *seq, size_t len, long max_errors,
    pattern_scan_best_t *b, pattern_scan_best_t *f) {
  assert(backwards != NULL && forwards != NULL);
  assert(seq != NULL || len == 0);
  assert(max_errors >= 0);
  if (backwards->len > PATTERN_SCAN_WORD ||
      forwards->len > PATTERN_SCAN_WORD) {
    /* (a column for each pattern: one after the other) */
    b->found = pattern_scan_best(backwards, seq, len, max_errors, true,
                                 &b->n_errors, &b->first, &b->last);
    f->found = pattern_scan_best(forwards, seq, len, max_errors, false,
                                 &f->n_errors, &f->first, &f->last);
    return;
  }
  /* the two scans do not depend on each other: in one loop, the
     processor computes them at the same time */
  const uint64_t *eq_b = backwards->peq, *eq_f = forwards->peq;
  /* (the last position of each pattern, as a mask: in this loop it
     takes less time than the shift of scan_step()) */
  const uint64_t top_b = (uint64_t)1 << (backwards->len - 1);
  const uint64_t top_f = (uint64_t)1 << (forwards->len - 1);
  uint64_t pv_b = ~(uint64_t)0, mv_b = 0, pv_f = ~(uint64_t)0, mv_f = 0;
  long score_b = (long)backwards->len, score_f = (long)forwards->len;
  long lowest_b = max_errors + 1, lowest_f = max_errors + 1;
  pattern_scan_best_t best_b = {0}, best_f = {0};
  for (size_t t = 0; t < len; t++) {
    size_t j = len - 1 - t;
    uint64_t up_b, down_b, up_f, down_f;
    word_step(eq_b[(unsigned char)seq[j]], &pv_b, &mv_b, &up_b, &down_b);
    word_step(eq_f[(unsigned char)seq[t]], &pv_f, &mv_f, &up_f, &down_f);
    score_b += (long)((up_b & top_b) != 0) - (long)((down_b & top_b) != 0);
    score_f += (long)((up_f & top_f) != 0) - (long)((down_f & top_f) != 0);
    note_score(score_b, (long)j, max_errors, true, &lowest_b, &best_b);
    note_score(score_f, (long)t, max_errors, false, &lowest_f, &best_f);
  }
  best_b.found = lowest_b <= max_errors;
  best_b.n_errors = lowest_b;
  best_f.found = lowest_f <= max_errors;
  best_f.n_errors = lowest_f;
  *b = best_b;
  *f = best_f;
}

/*
   Where the match of the pattern which ends at <end> of <seq> with the
   fewest errors (at most <max_errors>) begins; of equally good starts,
   the first (the longest match). The pattern is aligned backwards from
   <end>: col[i] is the distance of its last i bases to the bases read.
*/
static long match_start(pattern_scan_t *self, const char *seq, long end,
    long max_errors, long *n_errors) {
  const long m = (long)self->len;
  long *col = self->col;
  long most = m + max_errors < end + 1 ? m + max_errors : end + 1;
  long best = m, best_t = 0;   /* (no base of the sequence: m errors) */
  for (long i = 0; i <= m; i++)
    col[i] = i;
  for (long t = 1; t <= most; t++) {
    unsigned char c = (unsigned char)seq[end - t + 1];
    long diagonal = col[0];
    col[0] = t;
    for (long i = 1; i <= m; i++) {
      long d = diagonal + (equal_at(self, (size_t)(m - i), c) ? 0 : 1);
      if (col[i] + 1 < d)
        d = col[i] + 1;
      if (col[i - 1] + 1 < d)
        d = col[i - 1] + 1;
      diagonal = col[i];
      col[i] = d;
    }
    if (col[m] <= best) {
      best = col[m];
      best_t = t;
    }
  }
  *n_errors = best;
  return end - best_t + 1;
}

size_t pattern_scan_matches(pattern_scan_t *self, const char *seq,
    size_t len, long max_errors, pattern_match_t *matches) {
  assert(matches != NULL);
  if (max_errors >= (long)self->len)
    max_errors = (long)self->len - 1;
  pattern_scan_end_t ends[PATTERN_SCAN_MAX_MATCHES];
  size_t n = pattern_scan_run(self, seq, len, max_errors, ends,
                              PATTERN_SCAN_MAX_MATCHES);
  /* (after the scan: a long pattern uses the same column) */
  for (size_t e = 0; e < n; e++) {
    matches[e].end = ends[e].best;
    matches[e].start = match_start(self, seq, ends[e].best, max_errors,
                                   &matches[e].n_errors);
    assert(matches[e].n_errors == ends[e].n_errors);
  }
  return n;
}

/* --- unit tests --- */

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

/* fill <seq> with <len> random characters of <alphabet> */
static void random_chars(char *seq, size_t len, const char *alphabet,
    unsigned long *state) {
  size_t n = strlen(alphabet);
  for (size_t i = 0; i < len; i++) {
    *state = *state * 6364136223846793005UL + 1442695040888963407UL;
    seq[i] = alphabet[(*state >> 33) % n];
  }
  seq[len] = '\0';
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
  test_eq_long(ends[0].best, 15);
  test_eq_long(ends[0].n_errors, 0);
  test_eq_long(ends[1].best, 35);
  /* with one error allowed, the positions next to the end are
     acceptable too: one run for each match */
  test_eq_size(pattern_scan_run(p, seq, strlen(seq), 1, ends, 4), 2);
  test_eq_long(ends[0].first, 14);
  test_eq_long(ends[0].best, 15);
  test_eq_long(ends[0].last, 16);
  /* no room for the second run: not reported */
  test_eq_size(pattern_scan_run(p, seq, strlen(seq), 1, ends, 1), 1);
  test_eq_long(ends[0].best, 15);
  /* a mismatch, a deletion, an N in the sequence; nothing in a sequence
     without the pattern, or in an empty one */
  seq = "CCAGCAATAGCAGGCCCCAGCAAAGCAGGCCCCAGCANAAGCAGG";
  test_eq_size(pattern_scan_run(p, seq, strlen(seq), 1, ends, 4), 3);
  test_eq_long(ends[0].n_errors, 1);
  test_eq_long(ends[1].n_errors, 1);
  test_eq_long(ends[2].n_errors, 0);
  test_eq_size(pattern_scan_run(p, "CCCCCCCCCCCCCCCCCCCC", 20, 2, ends, 4),
               0);
  test_eq_size(pattern_scan_run(p, "", 0, 2, ends, 4), 0);
  /* the matches themselves: where they begin as well */
  pattern_match_t matches[PATTERN_SCAN_MAX_MATCHES];
  test_eq_size(pattern_scan_matches(p, seq, strlen(seq), 1, matches), 3);
  test_eq_long(matches[0].start, 2);
  test_eq_long(matches[0].end, 13);
  test_eq_long(matches[0].n_errors, 1);
  test_eq_long(matches[1].start, 18);
  test_eq_long(matches[1].end, 28);
  test_eq_long(matches[2].start, 33);
  test_eq_long(matches[2].end, 44);
  test_eq_long(matches[2].n_errors, 0);
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
      for (size_t i = 0; i < m; i++) {
        state = state * 6364136223846793005UL + 1442695040888963407UL;
        pattern[i] = "ACGTRN"[(state >> 33) % 6];
      }
      pattern[m] = '\0';
      size_t n = 0;
      while (n + m + 20 < sizeof(text)) {
        for (size_t i = 0; i < 20 + m; i++) {
          state = state * 6364136223846793005UL + 1442695040888963407UL;
          unsigned long r = state >> 33;
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
