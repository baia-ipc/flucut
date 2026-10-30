/* (c) 2026 Giorgio Gonnella */

/*
   Search of a primer in a read, using edlib (edit distance, the primer
   must be completely aligned, the ends of the read are free).

   Two points are specific to the primers searched here:

   - degenerate bases: the primers contain IUPAC codes (e.g. R = A or G,
     see primers.h), which match any of the bases they stand for; the
     same holds for the codes in the read (e.g. N);
   - the position of the conserved terminus in the read: the bases of a
     read which lie under a primer come from the primer, not from the
     template (see trimming.h), but only the terminus part of the primer
     is also template sequence. The trimming therefore needs to know
     where, in the read, the terminus of the primer begins. It is
     computed from the alignment path, so that insertions and deletions
     in the read move the boundary accordingly.

   If a partial alignment is allowed on one side and the complete primer
   is not found, the primer may overhang the read on that side, by up to
   half of its length; the overhanging part is not counted as errors, and
   the maximum number of errors is reduced proportionally to the part of
   the primer inside the read. A read can therefore start inside a
   primer, as reads of an amplicon often do.

   Among equally good locations, the longest one is preferred; the
   aligned region is then extended as long as the edit distance does not
   increase (i.e. primer bases unaligned at the ends of the alignment are
   counted as mismatches instead).
*/

#ifndef PRIMER_SEARCH_H
#define PRIMER_SEARCH_H

#include <stdbool.h>
#include <stddef.h>
#include "primers.h"

typedef struct primer_search_t primer_search_t;

/* result of a search */
typedef struct {
  /* start and end positions of the alignment on the read (0-based,
     inclusive); -1, -1 if not found */
  long start, end;
  /* the primer overhangs the read */
  bool partial;
  /* number of errors (if found) */
  long n_errors;
  /*
     Position in the read of the first base of the conserved terminus of
     the primer (the first base of the template-derived part, see
     primers.h). It is not clipped to the read: it is negative if the
     read begins inside the terminus, and at least the read length if the
     terminus is entirely after the end of the read. With a primer which
     is only a terminus (no tail), it is the start of the alignment.
  */
  long terminus_start;
} primer_hit_t;

/* Create a new searcher (it contains buffers: use one for each thread) */
primer_search_t* primer_search_new(void);

/*
   Search the primer <primer> in the read <read> (length <read_len>),
   with at most <max_errors> errors; <partial_left> and <partial_right>:
   whether the primer may overhang the read on that side.

   Returns true if found; the result is stored in <hit>.
*/
bool primer_search_run(primer_search_t *self, const primer_t *primer,
    const char *read, size_t read_len, long max_errors, bool partial_left,
    bool partial_right, primer_hit_t *hit);

/*
   The same, for a sequence which is not a primer of the built-in table
   (e.g. a tail, a barcode or a sequencing adapter): no terminus is
   computed (<terminus_start> is set to the end of the alignment plus
   one, i.e. the whole alignment is tail).
*/
bool primer_search_run_seq(primer_search_t *self, const char *seq,
    size_t seq_len, const char *read, size_t read_len, long max_errors,
    bool partial_left, bool partial_right, primer_hit_t *hit);

/* Delete a searcher (NULL is accepted) */
void primer_search_delete(primer_search_t *self);

/* Unit tests */
void primer_search_test(void);

#endif
