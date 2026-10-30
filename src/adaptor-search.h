/* (c) 2026 Giorgio Gonnella */

/*
   Search of an adaptor in a read, using edlib (edit distance, the
   adaptor must be completely aligned, the ends of the read are free).
   From Lumipore (github.com/baia-ipc/lumipore, src/adaptor-search.h),
   unchanged. An N of the adaptor matches any base (e.g. a barcode of
   unknown sequence, see adapters.h).

   If a partial alignment is allowed on one side and the complete adaptor
   is not found, the adaptor may overhang the read on that side, by up to
   half of its length; the overhanging part is not counted as errors, and
   the maximum number of errors is reduced proportionally to the part of
   the adaptor inside the read.

   Among equally good locations, the longest one is preferred; the aligned
   region is then extended as long as the edit distance does not increase
   (i.e. adaptor bases unaligned at the ends of the alignment are counted
   as mismatches instead).
*/

#ifndef ADAPTOR_SEARCH_H
#define ADAPTOR_SEARCH_H

#include <stdbool.h>
#include <stddef.h>

typedef struct adaptor_search_t adaptor_search_t;

/* result of a search */
typedef struct {
  /* start and end positions of the alignment on the read (0-based,
     inclusive); -1, -1 if not found */
  long start, end;
  /* the adaptor overhangs the read */
  bool partial;
  /* number of errors (if found) */
  long n_errors;
} adaptor_hit_t;

/* Create a new searcher (it contains buffers: use one for each thread) */
adaptor_search_t* adaptor_search_new(void);

/*
   Search the adaptor <adaptor> (length <adaptor_len>) in the read <read>
   (length <read_len>), with at most <max_errors> errors; <partial_left> and
   <partial_right>: whether the adaptor may overhang the read on that side.

   Returns true if found; the result is stored in <hit>.
*/
bool adaptor_search_run(adaptor_search_t *self, const char *adaptor,
    size_t adaptor_len, const char *read, size_t read_len, long max_errors,
    bool partial_left, bool partial_right, adaptor_hit_t *hit);

/* Delete a searcher (NULL is accepted) */
void adaptor_search_delete(adaptor_search_t *self);

/* Unit tests */
void adaptor_search_test(void);

#endif
