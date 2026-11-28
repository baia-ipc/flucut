/* (c) 2026 Giorgio Gonnella */

/*
   The universal primers of the influenza multisegment RT-PCR, built in:
   Uni12/Uni13 for influenza A and the IBV-GA2 cocktail for influenza B.

   The primers themselves - their sequences, the conserved termini they
   bind, the segments each one binds and the publication each is quoted
   from - are in metadata/primers.yaml, which the build turns into the
   tables used here. What is here is what the program does with them.

   The primers bind the conserved termini of the segments, so every
   amplicon begins and ends with primer-derived sequence: the bases
   under a primer reflect the primer and not the template, and are
   removed (see trimming.h). A read can be sequenced in either
   orientation, so both strands are searched (see primer-search.h).

   A primer which reaches past the universal terminus, as those of
   influenza B do, binds only the segments it was designed for, and its
   bases are the template only at an end of one of them: which segments
   those are is part of each primer (see primer_t.segments).

   Degenerate bases are written in the IUPAC code (R = A or G); they are
   matched against any of the bases they stand for (see iupac_bases() in
   iupac.h).
*/

#ifndef PRIMERS_H
#define PRIMERS_H

#include <stdbool.h>
#include <stddef.h>
#include "segments.h"
#include "virus.h"

/*
   The end of the segment at which a primer binds: PRIMER_END_LEFT is
   the terminus at the 3' end of the vRNA, which is the start of the
   plus strand of the amplicon (12 nt in influenza A, 9 in influenza
   B), PRIMER_END_RIGHT the terminus at the 5' end of the vRNA, the end
   of the plus strand (13 nt in influenza A, 10 in influenza B).
*/
typedef enum {
  PRIMER_END_LEFT,
  PRIMER_END_RIGHT,
} primer_end_t;

/*
   The segments whose terminus a primer binds, as a bit per segment (see
   segments.h): PRIMER_SEGMENTS_ALL for a primer which binds the same
   terminus in every segment, as the universal primers of influenza A
   do.

   A primer which reaches past the universal terminus into bases which
   differ between the segments binds only the segments it was designed
   for, and its 3' bases are the terminus of those segments alone. The
   bases it writes into a read are then the template only if the read
   comes from one of them: at an end of another segment they are the
   primer's own sequence, as at an internally primed end, and are
   removed (see amplicon.h).
*/
typedef unsigned primer_segments_t;
#define PRIMER_SEGMENTS_ALL ((primer_segments_t)0)
#define PRIMER_SEGMENT(S) ((primer_segments_t)1 << ((S) - 1))

/*
   A primer: the oligonucleotide as it is ordered, 5' -> 3', and the
   number of its trailing bases which are the conserved terminus of the
   segment (the rest, at the 5' end, is the tail added for the protocol).

   The tail and the terminus can overlap in the oligonucleotide: in the
   primers of Zhou et al. (2009) the last two bases of the published
   tail (ACGCGTGATCAG) are also the first two bases of the conserved
   terminus. The oligonucleotide is therefore stored whole, and
   <terminus_len> says where the template-derived part begins, instead of
   joining a tail and a terminus.

   <terminus_from> is the offset of the terminus in the
   oligonucleotide: in a primer as it is ordered the terminus is the
   suffix (<terminus_from> = <seq_len> - <terminus_len>), and in its
   reverse complement, searched at the other end of an amplicon, it is
   the prefix (<terminus_from> = 0, see primers_reverse()).
*/
typedef struct {
  const char *name;
  const char *seq;
  size_t seq_len;
  size_t terminus_from, terminus_len;
  primer_end_t end;
  /* the type of influenza the primer was designed for */
  virus_type_t virus;
  /* the segments whose terminus it binds (see above) */
  primer_segments_t segments;
} primer_t;

/*
   The conserved termini and the tails the primers are made of, and the
   built-in primers themselves, written by the build from
   metadata/primers.yaml (scripts/primers-table.awk). The sequences, the
   segments each primer binds and the publication each is quoted from
   are in that file; a primer or a scheme is added or changed there and
   not here.

   PRIMERS_TERMINUS_12 and PRIMERS_TERMINUS_13 are the conserved termini
   of influenza A, PRIMERS_TERMINUS_12_R the first with the degenerate
   fourth position; PRIMERS_TERMINUS_B_9 and PRIMERS_TERMINUS_B_10 (and
   PRIMERS_TERMINUS_B_10_T, its form in the HA, the NA and the NS) are
   those of influenza B.
*/
#include "primers-table.h"

/*
   The built-in primers, by name (case is ignored); NULL if there is no
   primer with that name.
*/
const primer_t* primers_builtin(const char *name);

/* The built-in primers, in order; <n> is set to their number */
const primer_t* primers_builtin_all(size_t *n);

/*
   Whether the terminus of the primer <p> is the terminus of segment
   <segment> (1 to 8, see segments.h), i.e. whether its bases are the
   template at an end of that segment. A primer of every segment
   (PRIMER_SEGMENTS_ALL) binds any of them; a segment-specific primer
   binds only its own, and an unknown segment (0) is not one of them,
   so that its bases are removed when the segment is not known.
*/
bool primers_binds_segment(const primer_t *p, int segment);

/*
   Write the reverse complement of <seq> (length <len>) to <dest>,
   complementing the degenerate codes as well (e.g. R -> Y); <dest> must
   have room for <len> + 1 characters. The result is 0-terminated. The
   sequence must be valid IUPAC (see iupac_is_seq() in iupac.h).
*/
void primers_reverse_complement(const char *seq, size_t len, char *dest);

/*
   The reverse complement of the primer <p>, as it appears at the other
   end of an amplicon: its sequence is written to <seq_buf> (which must
   have room for <p>->seq_len + 1 characters and must exist as long as
   <dest> is used), and its terminus becomes its prefix. The name is that
   of <p>.
*/
void primers_reverse(const primer_t *p, char *seq_buf, primer_t *dest);

/* Unit tests */
void primers_test(void);

#endif
