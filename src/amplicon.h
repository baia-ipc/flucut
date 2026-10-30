/* (c) 2026 Giorgio Gonnella */

/*
   Analysis of a read of the multisegment RT-PCR and computation of the
   part of it to keep.

   A read is searched for the primers of the scheme at both of its ends
   (see scheme.h). As a read can be sequenced in either orientation, the
   search is made on the read and, if the layout is not complete, on its
   reverse complement; the orientation with the better result is kept.
   The positions reported are always those of the original read, which is
   not reverse-complemented: a trimmer must not change the orientation of
   the reads it is given.

   The trimming
   ------------

   The universal primers bind the conserved termini of the segments, so
   the bases of a read which lie under a primer come from the primer and
   not from the template: a difference between the template and the
   primer at such a position is overwritten by the primer, and a variant
   there is an artefact. The part of the primer which is a tail added by
   the protocol is not template sequence at all.

   Because flucut knows which primers were used, it knows where, inside
   the primer, the tail ends and the conserved terminus begins (see
   primers.h). It can therefore offer a choice which a trimmer that only
   matches a sequence cannot:

   - TRIM_PRIMER (the default): all the primer-derived bases are removed,
     tail and conserved terminus together. No base of the output comes
     from a primer, so no variant of the consensus at the segment ends
     can come from one. The consensus is then shorter than the segment
     by the length of the termini at each end (12 and 13 bases), which
     matters for the programs which expect complete termini (e.g. the
     FLU-utr module of IRMA): the ends are either missing or have to be
     padded.
   - TRIM_TERMINUS: only the tail (and anything outside the primer) is
     removed, and the bases of the conserved terminus are kept. The
     consensus then has complete termini, but the kept bases are those of
     the primer, not those of the template.

   The removal is by position, not by difference: the whole span of the
   primer alignment is removed, whether or not the bases under it differ
   from the primer. A read which begins inside a primer is handled the
   same way, as the search allows a primer to overhang the read.
*/

#ifndef AMPLICON_H
#define AMPLICON_H

#include <stdbool.h>
#include <stddef.h>
#include "primer-search.h"
#include "scheme.h"

/* what is removed at the conserved ends (see above) */
typedef enum {
  TRIM_PRIMER,
  TRIM_TERMINUS,
} trim_policy_t;

/* the orientation in which the layout of the amplicon was found */
typedef enum {
  STRAND_PLUS,
  STRAND_MINUS,
  STRAND_NONE,     /* no primer was found in either orientation */
} strand_t;

/* the primer found at one end of a read */
typedef struct {
  bool found;
  size_t variant;      /* which primer of that end of the scheme */
  primer_hit_t hit;    /* positions in the analyzed orientation */
} amplicon_end_t;

typedef struct {
  strand_t strand;
  /* the ends of the amplicon: left and right in the orientation of the
     amplicon, i.e. of the analyzed sequence (for a read of the minus
     strand, left is the end of the original read) */
  amplicon_end_t left, right;
  /* the part of the ORIGINAL read to keep, 0-based and inclusive;
     <keep_to> < <keep_from> if nothing is left (e.g. a primer dimer) */
  long keep_from, keep_to;
} amplicon_t;

typedef struct amplicon_analyzer_t amplicon_analyzer_t;

/*
   Create an analyzer (it contains buffers: use one for each thread);
   <scheme> must exist as long as the analyzer is used.

   <max_errors>: the errors allowed in a primer; a number of bases, or,
   if smaller than 1, a fraction of the length of the primer.
*/
amplicon_analyzer_t* amplicon_analyzer_new(const scheme_t *scheme,
    double max_errors, trim_policy_t policy);

/*
   Analyze the read <seq> (length <len>); the result is valid until the
   next call on the same analyzer.
*/
const amplicon_t* amplicon_analyzer_run(amplicon_analyzer_t *self,
    const char *seq, size_t len);

/* Delete an analyzer (NULL is accepted) */
void amplicon_analyzer_delete(amplicon_analyzer_t *self);

/* Number of ends at which a primer was found (0, 1 or 2) */
static inline size_t amplicon_n_ends(const amplicon_t *self) {
  return (size_t)(self->left.found ? 1 : 0) +
         (size_t)(self->right.found ? 1 : 0);
}

/* Length of the part of the read which is kept */
static inline size_t amplicon_kept_length(const amplicon_t *self) {
  return self->keep_to < self->keep_from ? 0 :
         (size_t)(self->keep_to - self->keep_from + 1);
}

/* Unit tests */
void amplicon_test(void);

#endif
