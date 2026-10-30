/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <string.h>
#include "sequtil.h"
#include "test.h"

/* complement of each character (IUPAC codes, U to A, case kept); 0 for
   the other characters, which are their own complement (see complement()) */
static const unsigned char complement_table[256] = {
  ['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A', ['U'] = 'A',
  ['R'] = 'Y', ['Y'] = 'R', ['K'] = 'M', ['M'] = 'K',
  ['B'] = 'V', ['V'] = 'B', ['D'] = 'H', ['H'] = 'D',
  ['a'] = 't', ['c'] = 'g', ['g'] = 'c', ['t'] = 'a', ['u'] = 'a',
  ['r'] = 'y', ['y'] = 'r', ['k'] = 'm', ['m'] = 'k',
  ['b'] = 'v', ['v'] = 'b', ['d'] = 'h', ['h'] = 'd',
};

static inline char complement(char c) {
  unsigned char t = complement_table[(unsigned char)c];
  return t != 0 ? (char)t : c;
}

void reverse_complement(const char *seq, size_t len, char *dest) {
  assert(seq != NULL);
  assert(dest != NULL);
  for (size_t i = 0; i < len; i++)
    dest[i] = complement(seq[len - 1 - i]);
  dest[len] = '\0';
}

void py_slice(long len, long start, long stop, long *from, long *slice_len) {
  assert(len >= 0);
  assert(from != NULL);
  assert(slice_len != NULL);
  if (start < 0) {
    start += len;
    if (start < 0)
      start = 0;
  } else if (start > len)
    start = len;
  if (stop < 0) {
    stop += len;
    if (stop < 0)
      stop = 0;
  } else if (stop > len)
    stop = len;
  *from = start;
  *slice_len = stop > start ? stop - start : 0;
}

void sequtil_test(void) {
  char dest[16];
  reverse_complement("ACGTNacgtU", 10, dest);
  test_assert(strcmp(dest, "AacgtNACGT") == 0);
  reverse_complement("RYKMBVDHSW", 10, dest);
  test_assert(strcmp(dest, "WSDHBVKMRY") == 0);
  reverse_complement("", 0, dest);
  test_assert(strcmp(dest, "") == 0);
  long from, slice_len;
  py_slice(10, 2, 5, &from, &slice_len);
  test_eq(from, 2L, long, "%ld");
  test_eq(slice_len, 3L, long, "%ld");
  py_slice(10, 8, 20, &from, &slice_len);
  test_eq(slice_len, 2L, long, "%ld");
  py_slice(10, -3, 10, &from, &slice_len);
  test_eq(from, 7L, long, "%ld");
  test_eq(slice_len, 3L, long, "%ld");
  py_slice(10, -3, 2, &from, &slice_len);
  test_eq(slice_len, 0L, long, "%ld");
  py_slice(10, 5, 5, &from, &slice_len);
  test_eq(slice_len, 0L, long, "%ld");
  py_slice(10, -20, -8, &from, &slice_len);
  test_eq(from, 0L, long, "%ld");
  test_eq(slice_len, 2L, long, "%ld");
}
