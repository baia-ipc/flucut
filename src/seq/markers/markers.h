/* (c) 2026 Giorgio Gonnella */

/*
   The amino-acid markers read from the reads.

   A marker is a position of a protein whose residue is worth knowing
   before anything is assembled: the residues that decide whether a
   virus binds the receptors of birds or of humans, replicates in a
   mammalian host, or resists a drug. What is reported is not a call on
   the virus: it is what the reads carry at that position, with how many
   of them carry it, which also shows a residue present in a part of the
   reads only - a minority that a consensus would round away.

   What is in here
   ---------------

   A marker is listed only if it fulfils both of these:

   - its evidence is settled: the residue has a published phenotype
     which does not change with the naming of the clades, so that the
     table needs no upkeep;
   - it can be found locally: a short conserved peptide lies directly
     beside it, so that the position is reached from the sequence
     itself.

   Clades and genotypes are therefore not here. They are defined by
   nomenclatures which are revised, which no table compiled into a
   program can follow.

   How a marker is found
   ---------------------

   By the anchor beside it (anchor-site.h): the bases of the reference
   which encode the conserved peptide. The anchor is searched in each
   read, in both orientations and with errors allowed, and the window
   beside the match is read as the read carries it. Nothing is aligned
   and nothing is assembled, so a marker is available as the reads
   arrive, and a read with an indel elsewhere still gives it.

   Numbering
   ---------

   The position of a residue depends on what it is counted on: the
   haemagglutinin is numbered on its precursor, on the mature protein or
   on the structure of H3, and the neuraminidase of N1 and of N2 differ
   by one. No marker here uses a position in the sequence: the residue
   is the one beside the anchor, whatever its number in the sample, and
   the number in <positions> says only which published residue it is,
   together with the convention in <numbering>. A virus whose stalk is
   shorter, or whose signal peptide is one residue longer, therefore
   still gives the right residue.
*/

#ifndef MARKERS_H
#define MARKERS_H

#include <stdbool.h>
#include <stddef.h>
#include "anchor-site.h"
#include "fasta.h"
#include "logger.h"
#include "reference.h"

/* the bases of an anchor taken from the reference (eight codons) */
#define MARKER_ANCHOR_NT 24
#define MARKER_ANCHOR_AA (MARKER_ANCHOR_NT / 3)

/* the longest window read beside an anchor, and the longest value a
   marker reports (see marker_kind_t) */
#define MARKER_WINDOW_NT 276
#define MARKER_VALUE_MAX 20

/* the residues of a C-terminus reported (MARKER_TO_STOP) */
#define MARKER_TERMINAL_AA 4

/* how the value of a marker is read beside its anchor */
typedef enum {
  /* the residues of the window: the marker is those residues */
  MARKER_RESIDUES,
  /* the last residues before the first stop codon after the anchor,
     the stop included ("ESEV*"): the last residues of a protein */
  MARKER_TO_STOP,
  /* the number of codons from the end of the anchor to the first stop
     codon after it: how long a reading frame is */
  MARKER_FRAME_LENGTH,
  /* the bases between the first base of the anchor and the first base
     of the second anchor: how long the stretch between them is */
  MARKER_SPAN,
  MARKER_N_KINDS
} marker_kind_t;

/* one marker of the table */
typedef struct {
  const char *name;        /* in the reports: "PB2-627" */
  const char *protein;     /* "PB2" */
  size_t segment;          /* 1 to 8 */
  /* the peptide the anchor encodes, as a pattern of
     codons_find_peptide() ('.' matches any residue); MARKER_ANCHOR_NT
     bases from its first codon are the anchor */
  const char *peptide;
  /* MARKER_SPAN: the peptide of the second anchor, else NULL */
  const char *peptide2;
  /* where the window lies beside the anchor: a negative offset is
     bases before its first base, a positive one or zero bases after
     its last base. The two are counted from the two ends of the match
     in the read, so that an indel inside the anchor moves neither */
  long offset;
  size_t window;           /* bases read (0 for MARKER_SPAN) */
  marker_kind_t kind;
  const char *numbering;   /* what the positions are counted on */
  const char *positions;   /* the published positions: "627" */
  /* the residues usually carried, and those worth looking at, each
     one residue per position of the window ("": none); a value is
     flagged when it has one of <attention> at the position where
     <attention> has it */
  const char *usual;
  const char *attention;
  /* MARKER_FRAME_LENGTH: a frame of fewer codons than this deserves
     attention; MARKER_SPAN: a stretch shorter than the one of the
     reference by at least this many bases does (0: neither) */
  size_t attention_below;
  /* a value deserving attention is also a remark on the sample (see
     qc-flags.h). False where nearly every virus now circulating
     carries that residue, so that the remark would be on every
     sample: the value is then reported without one */
  bool flag;
  const char *meaning;     /* one line, for the reports */
} marker_def_t;

/*
   The markers and their number (MARKER_PB2_627, ..., MARKER_N_MARKERS),
   and, where MARKERS_TABLE_DEFINE is defined, the table itself: all of
   it written by the build from metadata/markers.yaml, so that a marker is
   added or changed in that file and not in the sources
*/
#include "markers-table.h"

/* The table, by marker (MARKER_PB2_627, ...) */
const marker_def_t* marker_def(size_t marker);

/*
   The anchors of one marker in one reference: the marker is read from
   the reads of a sample only if its reference has them
*/
typedef struct {
  bool found;              /* false: the other fields are not set */
  anchor_t anchor;
  anchor_t anchor2;        /* MARKER_SPAN */
  /* what the reference itself carries, which the reads are compared
     with, and the anchor translated, which the reads are compared with
     to see that it matched where it should */
  char value[MARKER_VALUE_MAX + 1];
  char anchor_aa[MARKER_ANCHOR_AA + 1];
} marker_reference_t;

/* the anchors of the markers in the references of a run, built once and
   shared */
typedef struct marker_anchors_t marker_anchors_t;

/*
   Build the anchors of the markers in <references> (NULL: none), whose
   segments <index> gives; NULL if no reference has any. An anchor is
   built only where the peptide beside the marker occurs exactly once in
   the reference and the bases of the anchor are found nowhere else in
   it within <max_error_rate> of its length as errors, which is how it
   is searched in a read: an anchor which another stretch of the segment
   could answer for would report the wrong position. A reference whose
   segment has markers but which gives none is reported at LOG_DEBUG and
   the marker is not read from it. The references must exist as long as
   the anchors are used.
*/
marker_anchors_t* marker_anchors_new(const fasta_t *references,
    const reference_index_t *index, double max_error_rate, Logger *logger);

/*
   The anchors to use for a sample, among the references <allowed>
   (NULL: all of them, see auto_settings_choose_references()): <out>,
   which must have room for MARKER_N_MARKERS pointers, is filled with
   the anchors of the first allowed reference of the segment of each
   marker which has them, or NULL where there is none. With the
   references of several subtypes, one reference per segment is allowed
   per sample, so these are the anchors of the subtype of that sample.
*/
void marker_anchors_chosen(const marker_anchors_t *self, const bool *allowed,
    const marker_reference_t **out);

/* How many (reference, marker) pairs have an anchor */
size_t marker_anchors_n(const marker_anchors_t *self);

void marker_anchors_delete(marker_anchors_t *self);

/*
   The value a window gives, into <out> (room for MARKER_VALUE_MAX + 1
   characters): the <window> bases of <nt> read as <kind> says. False,
   and <out> empty, if the window does not give one (no stop codon in it
   for MARKER_TO_STOP and MARKER_FRAME_LENGTH).
*/
bool marker_value_of(marker_kind_t kind, const char *nt, size_t window,
    char *out);

/* The value of a span of <bases>, as the number of bases (MARKER_SPAN) */
void marker_span_value(long bases, char *out);

/*
   Whether <value> deserves attention: a residue of <attention> at the
   position where <attention> has it; a reading frame of fewer codons
   than <attention_below>; or a stretch shorter than the one of the
   reference, <reference> ("": not known), by at least
   <attention_below> bases
*/
bool marker_value_attention(const marker_def_t *def, const char *value,
    const char *reference);

#endif
