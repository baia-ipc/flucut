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

   - TRIM_PRIMER: all the primer-derived bases are removed, tail and
     conserved terminus together. No base of the output comes from a
     primer, so no variant of the consensus can come from one, at the
     segment ends or at an internal priming site, and no reference is
     needed. The consensus is shorter than the segment by the termini
     (12 and 13 bases): those bases would have been the primer's, not
     measured sequence, and they are outside the coding sequences.
   - TRIM_TERMINUS: the conserved terminus is kept, but only at an end
     confidently classified as terminal (see below); the tail and
     anything outside the primer are removed. This gives complete segment
     ends while a primer base never reaches the consensus at an internal
     priming site. It needs the references (to classify the ends).

   The program chooses the mode from the references (options.c): with
   references, TRIM_TERMINUS (complete ends, and the internal priming
   located); without, TRIM_PRIMER. There is no option to override this.
   The kept conserved-terminus bases are those of the primer (with the
   degenerate Uni12 primer, possibly the base of the other variant);
   they are the segment's own conserved ends, the same in all segments
   and subtypes.

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

   A primer must be searched only near the ends of a read. A match deep
   inside a read is not the primer which started the amplicon, and
   trimming there would remove real template. The primer of an end is
   therefore only accepted within <window> bases of that end.

   The primer-derived bases of an internally primed read are at the end
   of the read, exactly as in a read of a real amplicon: the two cannot
   be told apart from the primer alone. TRIM_PRIMER removes the primer
   from both. TRIM_TERMINUS keeps the conserved terminus only where the
   classification against the references says the end is a segment
   terminus, and cuts the whole primer everywhere else, so the terminus
   of a genuine amplicon is kept while the footprint of an internally
   primed read is removed.

   Classification against the references
   -------------------------------------

   With the references of the segments (see reference.h), the bases
   which follow the primer of an end (the probe: for the right end, the
   reverse complement of the bases before the primer) are located in the
   segments. If they map at a terminus, the amplicon started there and
   the end is trimmed by the policy (with TRIM_TERMINUS, the conserved
   terminus is kept, so that the consensus has complete ends). If they
   map to an interior position, the amplicon was internally primed, and
   the whole primer of that end is removed: only that prefix, the rest of
   the read being genuine template. An end whose probe is not located (a
   subtype absent from the references, a probe too short) is reported as
   unknown and cut when unsure: its whole primer is removed, as at an
   internal end, so a primer base reaches the consensus only at an end
   known to be terminal.

   With TRIM_PRIMER (no references) the classification is not done; the
   whole primer is removed at every end.

   The matches found between the two end primers are reported as
   internal matches: they are the sites at which the primer can anneal.
   They are counted, not trimmed. The reads which contain several
   amplicons are split before the analysis (see chimera.h), so that the
   primers at their junctions are at the ends of the parts analyzed.

   A primer of the same end of the scheme can be found at both ends of a
   read, at one of them reverse-complemented: two amplicons joined
   without a primer between them, each read from its primer towards the
   junction. Only one of the two primers belongs to the layout found;
   the other is removed whole, with what lies outside it.

   Adapters
   --------

   At an end of the read at which no primer is found (a broken read, or
   one which begins too far inside the primer), the sequencing adapters
   of Oxford Nanopore are searched (see adapters.h): if one is found,
   it is removed with everything outside it. The part of the read
   between the adapters is then analyzed again: a primer whose outer
   bases are missing is now at an end of it, where it may overhang
   (with up to half of the primer missing), and is removed as well; the
   primer must then lie next to the adapter (within
   AMPLICON_ADAPTER_WINDOW bases), as a primer-like site within the
   usual window of the new end would otherwise be taken for it. The
   result of the second analysis is taken if it finds at least as many
   primers as the first.

   Placement
   ---------

   With references, a read of which no end was located (no primer at
   either end, or primers whose probes did not map) is placed on the
   references by a probe from inside its kept part (see placement.h).
   The placement is approximate and does not change the trimming: it is
   reported, so that the read is counted in the coverage of its segment.

   Ends located from inside
   ------------------------

   With references, an end of a read can be without a location: no
   primer at it (a broken molecule), or a primer whose probe did not map
   (the first bases after a primer are an untranslated region and the
   start of the coding sequence, where a sample can differ much from a
   reference). The span of such a read on the reference cannot be told
   from its length: a read with an internal deletion is much shorter
   than the distance between its ends. An end which was not located is
   therefore located by a probe from further inside the read, read from
   that end inwards; up to AMPLICON_INNER_PROBES probes are tried, each
   one probe length further in. If the other end has a position, the
   probe is searched in its reference, on the other strand (as the
   probes of the two ends of an amplicon are), by aligning it to the
   whole reference; otherwise by its k-mers, in all references. The position
   of the end is that of the probe, moved out by its distance from the
   end; it is approximate (it ignores the indels of the bases skipped).
   It does not change the origin of the end nor the trimming: it is
   reported, for the coverage (the span of the read, and whether it has
   an internal deletion).

   Internal deletions
   ------------------

   With references, a read whose two ends map to the same segment (an
   end located from inside included), much farther apart than the
   template between them is long (see junction_expected()), has an
   internal deletion (e.g. a defective interfering genome). Its junction
   is located by probes from inside the template (see junction.h). It
   does not change the trimming: it is reported, so that each flank of
   the read is counted where it lies in the coverage of its segment.
*/

#ifndef AMPLICON_H
#define AMPLICON_H

#include <stdbool.h>
#include <stddef.h>
#include "adapters.h"
#include "junction.h"
#include "placement.h"
#include "primer-search.h"
#include "reference.h"
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

/* the mode without references (the program derives the mode from the
   references: TRIM_TERMINUS with them, TRIM_PRIMER without, see above) */
#define AMPLICON_DEFAULT_POLICY TRIM_PRIMER

/* the fewest bases of a primer searched at an end of a read, after a
   flank, without its outer bases (see amplicon.c, missing_costs()) */
#define AMPLICON_MIN_TRUNCATED 8

/* default window: the primer of an end must begin (left) or end (right)
   within this many bases of that end of the read */
#define AMPLICON_DEFAULT_WINDOW 100

/* the greatest number of probes tried to locate an end from inside the
   read (see "Ends located from inside" above) */
#define AMPLICON_INNER_PROBES 4

/* the window of a primer next to an adapter (see "Adapters" above) */
#define AMPLICON_ADAPTER_WINDOW 5

/* a match of a primer inside a read (see above), in the positions of the
   original read */
typedef struct {
  bool at_left_end;    /* a primer of the left end of the scheme */
  size_t variant;      /* which primer of that end */
  long start, end;
  long n_errors;
} amplicon_match_t;

/* where the amplicon started, at one end (see "Classification against
   the references" above) */
typedef enum {
  ORIGIN_UNKNOWN,      /* no references, or the probe was not located */
  ORIGIN_TERMINAL,     /* at a terminus of a segment: a normal amplicon */
  ORIGIN_INTERNAL,     /* at an interior position: internal priming */
} end_origin_t;

/* the primer found at one end of a read */
typedef struct {
  bool found;
  size_t variant;      /* which primer of that end of the scheme */
  primer_hit_t hit;    /* positions in the analyzed orientation */
  end_origin_t origin;
  /* whether the conserved terminus of the primer found here is the
     terminus of the segment the probe was located in, so that its
     bases are the template and may be kept (TRIM_TERMINUS): the end is
     terminal and the primer binds that segment (see primers.h). False
     without references, as no end is then classified. */
  bool keep_terminus;
  reference_hit_t ref; /* the location of the probe (if located) */
  /* the position of the end from a probe further inside the read, if
     the end itself was not located (see "Ends located from inside"
     above); also for an end without a primer */
  reference_hit_t inner;
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
  /* the adapters found at the start and at the end of the ORIGINAL read
     (searched only at an end without a primer, see "Adapters" above) */
  adapter_hit_t adapter_start, adapter_end;
  /* where the kept part of the ORIGINAL read lies on the references, if
     no end was located (see "Placement" above) */
  placement_t placement;
  /* the junction of a read with an internal deletion, if located (see
     "Internal deletions" above) */
  junction_t junction;
} amplicon_t;

/* The position of an end of <a> on the references: that of its probe,
   or else that found from inside the read; NULL if neither */
static inline const reference_hit_t* amplicon_end_position(
    const amplicon_end_t *end) {
  return end->found && end->ref.mapped ? &end->ref :
         end->inner.mapped ? &end->inner : NULL;
}

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

/*
   Classify the ends of the reads against the references <refs> (see
   above); <refs> must exist as long as the analyzer is used, and belongs
   to the thread of the analyzer (it contains buffers). <terminal_window>:
   see reference_hit_is_terminal().
*/
void amplicon_analyzer_use_references(amplicon_analyzer_t *self,
    reference_set_t *refs, long terminal_window);

/*
   Search the adapters <adapters> at the ends of the reads at which no
   primer is found (see above); <adapters> must exist as long as the
   analyzer is used, and belongs to the thread of the analyzer.
*/
void amplicon_analyzer_use_adapters(amplicon_analyzer_t *self,
    adapters_t *adapters);

/* Classify the reads only against the references <allowed> (one value
   for each reference, see reference_set_restrict); NULL considers all.
   Without references it has no effect */
void amplicon_analyzer_restrict_references(amplicon_analyzer_t *self,
    const bool *allowed);

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
