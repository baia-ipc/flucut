/* (c) 2026 Giorgio Gonnella */

/*
   The universal primers of the influenza A multisegment RT-PCR
   (Uni12/Uni13), built in.

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

   Degenerate bases are written in the IUPAC code (R = A or G); they are
   matched against any of the bases they stand for (see
   primers_iupac_bases()).
*/

#ifndef PRIMERS_H
#define PRIMERS_H

#include <stdbool.h>
#include <stddef.h>

/*
   The end of the segment at which a primer binds: PRIMER_END_12 is the
   conserved 12 nt terminus (the 3' end of the vRNA, the start of the
   plus strand of the amplicon), PRIMER_END_13 the conserved 13 nt
   terminus (the 5' end of the vRNA, the end of the plus strand).
*/
typedef enum {
  PRIMER_END_12,
  PRIMER_END_13,
} primer_end_t;

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
*/
typedef struct {
  const char *name;
  const char *seq;
  size_t seq_len;
  size_t terminus_len;
  primer_end_t end;
} primer_t;

/* The conserved termini themselves (without any tail) */
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
*/
const primer_t* primers_builtin(const char *name);

/* The built-in primers, in order; <n> is set to their number */
const primer_t* primers_builtin_all(size_t *n);

/*
   The bases a IUPAC code stands for, as a 0-terminated string of
   A, C, G and T (e.g. "AG" for 'R', "A" for 'A'); NULL if <code> is not
   a IUPAC code. Case is ignored.
*/
const char* primers_iupac_bases(char code);

/* Whether a sequence contains only IUPAC codes (and is not empty) */
bool primers_is_iupac_seq(const char *seq, size_t len);

/*
   Append to <dest> the reverse complement of <seq> (length <len>),
   complementing the degenerate codes as well (e.g. R -> Y); <dest> must
   have room for <len> characters. The sequence must be valid (see
   primers_is_iupac_seq()).
*/
void primers_reverse_complement(const char *seq, size_t len, char *dest);

/* Unit tests */
void primers_test(void);

#endif
