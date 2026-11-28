/* (c) 2026 Giorgio Gonnella */

/*
   The universal primers of the influenza multisegment RT-PCR, built in:
   Uni12/Uni13 for influenza A and the IBV-GA2 cocktail for influenza B.

   All influenza A genomic segments end in the same conserved sequences:
   the 12 nucleotides of the 3' end and the 13 nucleotides of the 5' end
   of the vRNA are (almost) identical in all eight segments and in all
   subtypes. The universal primers of the multisegment RT-PCR bind these
   termini, so that every amplicon begins and ends with primer-derived
   sequence (see trimming.h: the bases under a primer reflect the primer,
   not the template).

   The sequences below are quoted from the publications cited for each of
   them; no sequence is derived or guessed. A sequence which cannot be
   verified in a publication or protocol is not added.

   In a double-stranded amplicon of a segment, the plus strand (the
   orientation of the mRNA) is

     [tail] AGCAAAAGCAGG ... segment ... CCTTGTTTCTACT [reverse tail]

   i.e. the Uni12 primer at the left end and the reverse complement of
   the Uni13 primer at the right end. A read can be sequenced in either
   orientation, so both strands are searched (see primer-search.h).

   The segments of influenza B end in conserved sequences as well, but
   shorter ones: 9 nucleotides at the start of the plus strand and 10 at
   its end. Nine and ten bases are too few for a specific primer, so the
   universal primers of influenza B reach past them, into the
   untranslated bases next to the terminus, which are conserved in a
   pattern shared by a few segments (Zhou et al. 2014, below). One
   primer therefore does not bind every segment, and its bases are the
   template only at an end of a segment it was designed for: which
   segments that is, is part of each primer (see primer_t.segments).

   The 10 nt terminus of influenza B has two forms, which differ at
   their sixth base: the HA, the NA and the NS segments have one, the
   other five segments the other. Both are built in, in the primers
   which use them.

   Degenerate bases are written in the IUPAC code (R = A or G); they are
   matched against any of the bases they stand for (see
   iupac_bases() in iupac.h).
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

/* The conserved termini of influenza A (without any tail) */
#define PRIMERS_TERMINUS_12 "AGCAAAAGCAGG"
#define PRIMERS_TERMINUS_13 "AGTAGAAACAAGG"
/*
   The same 12 nt terminus with the degenerate fourth position (A or G),
   which accommodates the variation of the segment termini
*/
#define PRIMERS_TERMINUS_12_R "AGCRAAAGCAGG"

/* The 5' tail shared by the primers of Zhou et al. (2009) */
#define PRIMERS_MBT_TAIL "ACGCGTGATCAG"

/*
   The conserved termini of influenza B: 9 bases at the start of the
   plus strand, 10 at its end, the last in the two forms which differ
   at the sixth base (PRIMERS_TERMINUS_B_10_T in the HA, the NA and the
   NS, PRIMERS_TERMINUS_B_10 in the other five segments). They are the
   first bases of the primers of Zhou et al. (2014) after the tail (see
   below).
*/
#define PRIMERS_TERMINUS_B_9 "AGCAGAAGC"
#define PRIMERS_TERMINUS_B_10 "AGTAGAAACA"
#define PRIMERS_TERMINUS_B_10_T "AGTAGTAACA"

/*
   The 5' tails of the IBV-GA2 primers of Zhou et al. (2014): six bases
   on the primers of the left end, ten on those of the right end. They
   are the promoter and terminator sequences of the reverse-genetics
   plasmid, not influenza sequence, and are always removed.
*/
#define PRIMERS_IBV_TAIL_LEFT "GGGGGG"
#define PRIMERS_IBV_TAIL_RIGHT "CCGGGTTATT"

/*
   The built-in primers, by name (case is ignored); NULL if there is no
   primer with that name.

   The names and the sequences:

   - "uni12"       AGCAAAAGCAGG
   - "uni13"       AGTAGAAACAAGG
     The conserved termini, used as primers without a tail. They are
     commonly called Uni12 and Uni13 after the universal primer set of
     Hoffmann et al. 2001, Arch Virol 146:2275-2289,
     https://doi.org/10.1007/s007050170002; the sequences used here are
     those quoted as Uni12 and Uni13 by Zhou et al. 2009 (see below),
     which are also the 3' ends of mbtuni-12 and mbtuni-13.
   - "uni12r"      AGCRAAAGCAGG
     The 12 nt terminus with the degenerate fourth position.
   - "mbtuni-12"   ACGCGTGATCAGCAAAAGCAGG
   - "mbtuni-13"   ACGCGTGATCAGTAGAAACAAGG
     The tailed primers of the single-reaction multisegment RT-PCR of
     Zhou et al. 2009, J Virol 83:10309-10313,
     https://doi.org/10.1128/JVI.01109-09. The common 5' tail
     ACGCGTGATCAG carries an MluI site (ACGCGT) used for cloning in that
     protocol.
   - "mbtuni-12.4r" ACGCGTGATCAGCRAAAGCAGG
     The tailed Uni12 primer with the degenerate fourth position, used
     in place of (or mixed with) mbtuni-12 by Licheri et al. 2025,
     Genome Med 17:103, https://doi.org/10.1186/s13073-025-01512-x
     (there MBTuni-12.4R; mixed with mbtuni-12 in a 1:4 ratio in the
     reverse transcription step).

   The primers of influenza B, with the segments each one binds:

   - "b-pbs-unif", "b-pbs-unir"     PB1 and PB2
   - "b-pa-unif", "b-pa-unir"       PA
   - "b-hana-unif", "b-hana-unir"   HA and NA
   - "b-np-unif", "b-np-unir"       NP
   - "b-m-uni3f", "b-mg-uni3f"      M (two forms of the primer of the
                                    left end, both in the cocktail)
   - "b-m-uni3r"                    M
   - "b-ns-uni3f", "b-ns-uni3r"     NS

     The universal IBV-GA2 primer cocktail of Zhou et al. 2014, J Clin
     Microbiol 52:1330-1337, https://doi.org/10.1128/JCM.03265-13
     (Table 1; the names are those of the table, in lower case). Six
     pairs amplify the eight segments: each primer is a 5' tail
     (PRIMERS_IBV_TAIL_LEFT or PRIMERS_IBV_TAIL_RIGHT, not influenza
     sequence) followed by the conserved terminus of its end and by
     five or eight further bases, which are the untranslated bases next
     to the terminus in the segments the primer binds.

     The cocktail mixes the primers in the proportions of Table 1; the
     proportions do not matter to flucut, which searches every primer
     of the scheme in every read.
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
