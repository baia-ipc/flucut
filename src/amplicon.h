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

   Internal priming
   ----------------

   A universal primer is short, and the untailed Uni12 primer is only
   the twelve bases of the conserved terminus, with a degenerate
   position. A sequence that similar occurs by chance inside the
   segments: the primer can anneal there and start an amplicon in the
   middle of a segment. The reads of such an amplicon begin with
   primer-derived bases which are not the template at that position, so
   they carry the sequence of the primer into the consensus, as a
   variant which looks real and is not.

   Two things follow for the trimming:

   - the primer-derived bases of an internally primed read are at the
     end of the read, exactly as in a read of a real amplicon, so
     TRIM_PRIMER removes them from both and the wrong bases never reach
     the consensus. TRIM_TERMINUS, which keeps the bases of the
     terminus, keeps the wrong bases: it must not be used when internal
     priming is possible;
   - a primer must be searched only near the ends of a read. A match
     deep inside a read is not the primer which started the amplicon,
     and trimming there would remove real template. The primer of an
     end is therefore only accepted within <window> bases of that end.

   The matches found between the two end primers are reported as
   internal matches: they are the sites at which the primer can anneal,
   and the junctions of chimeras and of reads containing several
   amplicons. They are counted, not trimmed.
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

/* the greatest number of internal matches kept for one read */
#define AMPLICON_MAX_INTERNAL 8

/* default window: the primer of an end must begin (left) or end (right)
   within this many bases of that end of the read */
#define AMPLICON_DEFAULT_WINDOW 100

/* a match of a primer inside a read (see above), in the positions of the
   original read */
typedef struct {
  bool at_left_end;    /* a primer of the left end of the scheme */
  size_t variant;      /* which primer of that end */
  long start, end;
  long n_errors;
} amplicon_match_t;

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
  /* the matches of a primer between the two ends (see above) */
  amplicon_match_t internal[AMPLICON_MAX_INTERNAL];
  size_t n_internal;
} amplicon_t;

typedef struct amplicon_analyzer_t amplicon_analyzer_t;

/*
   Create an analyzer (it contains buffers: use one for each thread);
   <scheme> must exist as long as the analyzer is used.

   <max_errors>: the errors allowed in a primer; a number of bases, or,
   if smaller than 1, a fraction of the length of the primer.

   <window>: the primer of an end is accepted only if it begins (left
   end) or ends (right end) within this many bases of that end of the
   read (see above; AMPLICON_DEFAULT_WINDOW).
*/
amplicon_analyzer_t* amplicon_analyzer_new(const scheme_t *scheme,
    double max_errors, trim_policy_t policy, long window);

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
