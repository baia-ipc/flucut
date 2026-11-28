/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "iupac.h"
#include "pattern-scan.h"

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
  const iupac_equal_chars_t *equal_chars = iupac_equal_chars();
  if (word)
    memset(self->peq, 0, sizeof(self->peq));
  else
    memset(self->equal, 0, 256 * len);
  for (size_t i = 0; i < len; i++) {
    const iupac_equal_chars_t *a = equal_chars + (unsigned char)pattern[i];
    for (int j = 0; j < a->n; j++) {
      unsigned char c = a->chars[j];
      if (word)
        self->peq[c] |= (uint64_t)1 << i;
      else
        self->equal[(size_t)c * len + i] = 1;
    }
  }
}

pattern_scan_t* pattern_scan_new(const char *pattern, size_t len) {
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

/* note the distance <score> at the position <j> of a scan for the
   best end (see pattern_scan_best()); <lowest>: the lowest so far */
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

bool pattern_scan_best(pattern_scan_t *self, const char *seq, size_t len,
    long max_errors, bool backwards, long *n_errors, long *first,
    long *last) {
  assert(self != NULL);
  assert(seq != NULL || len == 0);
  assert(max_errors >= 0);
  scan_state_t st;
  scan_start(self, &st);
  long lowest = max_errors + 1;
  pattern_scan_best_t best = {0};
  for (size_t t = 0; t < len; t++) {
    size_t j = backwards ? len - 1 - t : t;
    long score = scan_step(self, &st, (unsigned char)seq[j]);
    note_score(score, (long)j, max_errors, backwards, &lowest, &best);
  }
  if (lowest > max_errors)
    return false;
  *n_errors = lowest;
  *first = best.first;
  *last = best.last;
  return true;
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
