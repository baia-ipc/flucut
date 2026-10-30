/* (c) 2026 Giorgio Gonnella */

/*
   Bounded edit distance of short sequences (e.g. indices, barcodes, UMIs,
   compared with the sequences of a list): dynamic programming in the
   diagonal band of width 2k + 1 around the main diagonal, stopping as soon
   as a row exceeds the maximum distance k; for sequences of at most
   SHORT_EDIT_DISTANCE_MAX bases (for longer sequences, use e.g. edlib).

   The functions are static inline: the comparison of a sequence with all
   the keys of a list is often in the innermost loop of a program, and the
   code is then compiled in place, as if written in the caller.
*/

#ifndef EDIT_DISTANCE_H
#define EDIT_DISTANCE_H

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>

/* maximum length of the sequences compared by the functions */
#define SHORT_EDIT_DISTANCE_MAX 64

/* internal (used by the functions below; not part of the API, the
   trailing underscore marks it, as the leading one is reserved): the rows
   of the dynamic programming (a: rows, b: columns) in <row> (at
   least blen + 1 values), with the values over k set to k + 1; with
   <prefix>, b can be longer than the end of the alignment (the band can
   end before the end of b); false if a row exceeds k (or, with <prefix>,
   the band is past the end of b) */
static inline bool short_edit_rows_(const char *a, size_t alen,
    const char *b, size_t blen, long k, bool prefix, long *row) {
  const long inf = k + 1;
  for (size_t j = 0; j <= blen; j++)
    row[j] = (long)j <= k ? (long)j : inf;
  for (size_t i = 1; i <= alen; i++) {
    size_t lo = (long)i - k > 1 ? (size_t)((long)i - k) : 1;
    size_t hi = i + (size_t)k < blen ? i + (size_t)k : blen;
    if (prefix && lo > blen)
      return false;
    long diag = row[lo - 1];
    row[lo - 1] = lo == 1 ? ((long)i <= k ? (long)i : inf) : inf;
    long row_min = row[lo - 1];
    for (size_t j = lo; j <= hi; j++) {
      long up = row[j];
      long best = diag + (a[i - 1] != b[j - 1]);
      if (up + 1 < best)
        best = up + 1;
      if (row[j - 1] + 1 < best)
        best = row[j - 1] + 1;
      if (best > inf)
        best = inf;
      diag = up;
      row[j] = best;
      if (best < row_min)
        row_min = best;
    }
    if (hi < blen)
      row[hi + 1] = inf;
    if (row_min > k)
      return false;
  }
  return true;
}

/*
   Edit distance of a and b (lengths at most SHORT_EDIT_DISTANCE_MAX), if
   at most k (k >= 0), otherwise k + 1
*/
static inline long short_edit_distance(const char *a, size_t alen,
    const char *b, size_t blen, long k) {
  /* the row on the stack */
  assert(alen <= SHORT_EDIT_DISTANCE_MAX && blen <= SHORT_EDIT_DISTANCE_MAX);
  long diff = (long)alen - (long)blen;
  if (diff > k || -diff > k)
    return k + 1;
  long row[SHORT_EDIT_DISTANCE_MAX + 1];
  if (!short_edit_rows_(a, alen, b, blen, k, false, row))
    return k + 1;
  return row[blen] <= k ? row[blen] : k + 1;
}

/*
   Distance of the best alignment of the whole <query> with a prefix of
   <target> (lengths at most SHORT_EDIT_DISTANCE_MAX; target not empty),
   if at most k (k >= 0), otherwise -1; <end> is set to the length of the
   prefix (the shortest among those with the best distance, as edlib in
   mode SHW)
*/
static inline long short_prefix_distance(const char *query, size_t qlen,
    const char *target, size_t tlen, long k, size_t *end) {
  /* the row on the stack */
  assert(qlen <= SHORT_EDIT_DISTANCE_MAX && tlen <= SHORT_EDIT_DISTANCE_MAX);
  long row[SHORT_EDIT_DISTANCE_MAX + 1];
  if (!short_edit_rows_(query, qlen, target, tlen, k, true, row))
    return -1;
  long best = k + 1;
  size_t from = (long)qlen - k > 0 ? (size_t)((long)qlen - k) : 0;
  for (size_t j = from; j <= tlen; j++)   /* the band of the last row */
    if (row[j] < best) {
      best = row[j];
      *end = j;
    }
  return best <= k ? best : -1;
}

/* Unit tests */
void edit_distance_test(void);

#endif
