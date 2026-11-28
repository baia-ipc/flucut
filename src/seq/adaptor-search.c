/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "edlib.h"
#include "adaptor-search.h"
#include "alloc.h"
#include "numutils.h"
#include "test.h"

#define PAD_CHAR 'N'

/* the padding character matches any base (also in the read) */
static const EdlibEqualityPair pad_equalities[] = {
  {PAD_CHAR, 'A'}, {PAD_CHAR, 'C'}, {PAD_CHAR, 'G'}, {PAD_CHAR, 'T'}};
#define N_PAD_EQUALITIES \
  ((int)(sizeof(pad_equalities) / sizeof(pad_equalities[0])))

struct adaptor_search_t {
  char *padded;
  size_t padded_cap;
};

adaptor_search_t* adaptor_search_new(void) {
  adaptor_search_t *self;
  XMALLOC(self, sizeof(adaptor_search_t));
  self->padded = NULL;
  self->padded_cap = 0;
  return self;
}

void adaptor_search_delete(adaptor_search_t *self) {
  if (self == NULL)
    return;
  free(self->padded);
  free(self);
}

/* edit distance between the adaptor and target[start..end] */
static long distance(const char *adaptor, size_t adaptor_len,
    const char *target, long start, long end) {
  EdlibAlignResult result = edlibAlign(adaptor, (int)adaptor_len,
      target + start, (int)(end - start + 1),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE,
                          pad_equalities, N_PAD_EQUALITIES));
  assert(result.status == EDLIB_STATUS_OK);
  long d = result.editDistance;
  edlibFreeAlignResult(result);
  return d;
}

/*
   Search the adaptor in <target> (the read, of length <read_len>, padded
   by <pad_left> characters on the left and <target_len> - <read_len> -
   <pad_left> on the right); returns false if not found; otherwise the
   positions in the read (clipped to the read), whether the adaptor
   overhangs the read, the number of errors and the number of adaptor
   positions inside the read (in <in_read>, unless NULL)
*/
static bool search(const char *adaptor, size_t adaptor_len,
    const char *target, size_t target_len, long pad_left, long read_len,
    long max_errors, adaptor_hit_t *hit, long *in_read) {
  if (target_len == 0)
    return false;
  EdlibAlignResult result = edlibAlign(adaptor, (int)adaptor_len,
      target, (int)target_len,
      edlibNewAlignConfig((int)max_errors, EDLIB_MODE_HW, EDLIB_TASK_LOC,
                          pad_equalities, N_PAD_EQUALITIES));
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
         distance(adaptor, adaptor_len, target, start - 1, end) <= n_errors)
    start--;
  while (end < (long)target_len - 1 &&
         distance(adaptor, adaptor_len, target, start, end + 1) <= n_errors)
    end++;
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
  return true;
}

static void not_found(adaptor_hit_t *hit) {
  hit->start = -1;
  hit->end = -1;
  hit->partial = false;
  hit->n_errors = -1;
}

bool adaptor_search_run(adaptor_search_t *self, const char *adaptor,
    size_t adaptor_len, const char *read, size_t read_len, long max_errors,
    bool partial_left, bool partial_right, adaptor_hit_t *hit) {
  assert(self != NULL);
  assert(adaptor != NULL);
  assert(adaptor_len > 0);
  assert(read != NULL);
  assert(hit != NULL);
  /* (unpadded: the whole alignment is in the read, in_read not needed) */
  if (search(adaptor, adaptor_len, read, read_len, 0, (long)read_len,
             max_errors, hit, NULL))
    return true;
  if (!partial_left && !partial_right) {
    not_found(hit);
    return false;
  }
  size_t max_overhang = adaptor_len / 2;
  size_t pad_left = partial_left ? max_overhang : 0;
  size_t pad_right = partial_right ? max_overhang : 0;
  size_t padded_len = pad_left + read_len + pad_right;
  XCHECKSPACEINC(self->padded, &self->padded_cap, padded_len + 1, 1024);
  memset(self->padded, PAD_CHAR, pad_left);
  memcpy(self->padded + pad_left, read, read_len);
  memset(self->padded + pad_left + read_len, PAD_CHAR, pad_right);
  self->padded[padded_len] = '\0';
  long in_read;
  if (!search(adaptor, adaptor_len, self->padded, padded_len, (long)pad_left,
              (long)read_len, max_errors, hit, &in_read) ||
      hit->n_errors > floordiv(max_errors * in_read, (long)adaptor_len)) {
    not_found(hit);
    return false;
  }
  return true;
}

/* run a search and compare the result with the expected values
   (start == -1: not found) */
static void check(adaptor_search_t *s, const char *adaptor, const char *read,
    long max_errors, bool left, bool right, long start, long end,
    bool partial, long n_errors) {
  adaptor_hit_t hit;
  bool found = adaptor_search_run(s, adaptor, strlen(adaptor), read,
                                  strlen(read), max_errors, left, right, &hit);
  test_assert(found == (start >= 0));
  test_eq(hit.start, start, long, "%ld");
  test_eq(hit.end, end, long, "%ld");
  if (found) {
    test_assert(hit.partial == partial);
    test_eq(hit.n_errors, n_errors, long, "%ld");
  }
}

void adaptor_search_test(void) {
  adaptor_search_t *s = adaptor_search_new();
  check(s, "AGCTTA", "GGAGCTTAGG", 1, false, false, 2, 7, false, 0);
  check(s, "AGCTTA", "GGAGGTTAGG", 1, false, false, 2, 7, false, 1);
  check(s, "ACGTACGT", "", 2, false, false, -1, -1, false, 0);
  check(s, "ACGTACGT", "", 2, true, false, -1, -1, false, 0);
  check(s, "ACGTACGT", "TACGT", 2, true, false, 0, 4, true, 0);
  check(s, "ACGTACGT", "GGACGTAC", 2, false, true, 2, 7, false, 2);
  check(s, "ACGTACGT", "GGGGGGACGTA", 2, false, true, 6, 10, true, 0);
  check(s, "ACGTACGT", "GGGGGGACGTA", 2, false, false, -1, -1, false, 0);
  adaptor_search_delete(s);
}
