/* (c) 2026 Giorgio Gonnella */

/*
   Sequence utilities: reverse complement and substrings with the
   semantics of the Python slices (useful for reproducing exactly the
   output of Python programs).
*/

#ifndef SEQUTIL_H
#define SEQUTIL_H

#include <stddef.h>

/*
   Write the reverse complement of <seq> (length <len>) to <dest>
   (at least <len> + 1 bytes; the result is 0-terminated).

   IUPAC codes are complemented, U is complemented to A, the case is kept
   and other characters are unchanged.
*/
void reverse_complement(const char *seq, size_t len, char *dest);

/*
   Compute the bounds of the Python slice s[<start>:<stop>] of a string
   of length <len> (negative values count from the end; values out of the
   range are clipped).

   Output parameters: <from>, the start position of the slice, and
   <slice_len>, its length (possibly 0).
*/
void py_slice(long len, long start, long stop, long *from, long *slice_len);

/* Unit tests */
void sequtil_test(void);

#endif
