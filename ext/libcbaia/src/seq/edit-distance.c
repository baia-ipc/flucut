/* (c) 2026 Giorgio Gonnella */

#include <stdlib.h>
#include "edit-distance.h"
#include "test.h"

/* --- unit tests --- */

#define MIN3(A, B, C) ((A) < (B) ? ((A) < (C) ? (A) : (C)) : \
                                   ((B) < (C) ? (B) : (C)))

/* full dynamic programming: edit distance (prefix false) or best
   alignment of a with a prefix of b, shortest prefix (prefix true) */
static long full_distance(const char *a, size_t alen, const char *b,
    size_t blen, bool prefix, size_t *end) {
  long d[SHORT_EDIT_DISTANCE_MAX + 1][SHORT_EDIT_DISTANCE_MAX + 1];
  for (size_t i = 0; i <= alen; i++)
    d[i][0] = (long)i;
  for (size_t j = 0; j <= blen; j++)
    d[0][j] = (long)j;
  for (size_t i = 1; i <= alen; i++)
    for (size_t j = 1; j <= blen; j++)
      d[i][j] = MIN3(d[i - 1][j - 1] + (a[i - 1] != b[j - 1]),
                     d[i - 1][j] + 1, d[i][j - 1] + 1);
  if (!prefix)
    return d[alen][blen];
  long best = d[alen][0];
  *end = 0;
  for (size_t j = 1; j <= blen; j++)
    if (d[alen][j] < best) {
      best = d[alen][j];
      *end = j;
    }
  return best;
}

static void random_seq(char *s, size_t len) {
  for (size_t i = 0; i < len; i++)
    s[i] = "ACGT"[rand() % 4];
}

void edit_distance_test(void) {
  test_eq(short_edit_distance("ACGT", 4, "ACGT", 4, 2), 0L, long, "%ld");
  test_eq(short_edit_distance("ACGT", 4, "AGGT", 4, 2), 1L, long, "%ld");
  test_eq(short_edit_distance("ACGT", 4, "ACT", 3, 2), 1L, long, "%ld");
  test_eq(short_edit_distance("ACGT", 4, "TTTT", 4, 2), 3L, long, "%ld");
  test_eq(short_edit_distance("ACGTAC", 6, "AC", 2, 2), 3L, long, "%ld");
  test_eq(short_edit_distance("ACGT", 4, "ACGT", 4, 0), 0L, long, "%ld");
  size_t end;
  test_eq(short_prefix_distance("ACGT", 4, "ACGTTTT", 7, 1, &end), 0L, long,
          "%ld");
  test_eq_size(end, 4);
  test_eq(short_prefix_distance("ACGT", 4, "AGTTTT", 6, 1, &end), 1L, long,
          "%ld");
  test_eq_size(end, 3);
  test_eq(short_prefix_distance("ACGT", 4, "TTTTTT", 6, 1, &end), -1L, long,
          "%ld");
  /* random sequences, compared with the full dynamic programming */
  srand(42);
  char a[SHORT_EDIT_DISTANCE_MAX], b[SHORT_EDIT_DISTANCE_MAX];
  size_t n_different = 0;
  for (size_t t = 0; t < 20000; t++) {
    size_t alen = (size_t)(rand() % 13), blen = (size_t)(rand() % 13);
    long k = rand() % 4;
    random_seq(a, alen);
    random_seq(b, blen);
    if (t % 2)  /* similar sequences */
      for (size_t i = 0; i < alen && i < blen; i++)
        if (rand() % 4)
          b[i] = a[i];
    long full = full_distance(a, alen, b, blen, false, NULL);
    long d = short_edit_distance(a, alen, b, blen, k);
    if (d != (full <= k ? full : k + 1))
      n_different++;
    if (blen == 0)
      continue;
    size_t full_end, short_end = 0;
    full = full_distance(a, alen, b, blen, true, &full_end);
    d = short_prefix_distance(a, alen, b, blen, k, &short_end);
    if (d != (full <= k ? full : -1) || (d >= 0 && short_end != full_end))
      n_different++;
  }
  test_eq_size(n_different, 0);
}
