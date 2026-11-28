/* (c) 2026 Giorgio Gonnella */

/*
   The cleavage site of the haemagglutinin, and where to find it in a
   read.

   The site is the last residues of HA1, where the protein is cut in
   two. In the highly pathogenic avian H5 and H7 viruses it is a run of
   basic residues (e.g. PQRERRRKR), encoded by an A/G-rich stretch; in
   the low pathogenic ones a single one (e.g. PQRETR). It is the residue
   that decides the pathotype, so a consensus which is wrong there
   reports the wrong virus.

   That stretch is exactly the kind of low-complexity sequence the reads
   get wrong, by an insertion or a deletion, and an assembler called
   from few reads can carry such an error into the consensus: the
   reading frames of HA1 and of HA2 then differ and the motif is
   nonsense, or, with three bases too many, the run of basic residues
   is silently one residue longer.

   What the reads carry there is counted in cleavage-counts.h, which
   also says how far it can be trusted. Nothing here is assembled.

   Where it is read from
   ---------------------

   HA2 begins directly after the cleavage site with the fusion peptide
   (GLFGAIAG and a few variants, see CLEAVAGE_FUSION_PEPTIDE), which
   every haemagglutinin has and which lies outside the low-complexity
   stretch. Those bases of the HA reference are the anchor
   (anchor-site.h); the cleavage site is the CLEAVAGE_SITE_NT bases of
   the read immediately before it, read back in whole codons from the
   first base of HA2, so the frame is set by the anchor and no indel
   elsewhere in the read can shift it.
*/

#ifndef CLEAVAGE_SITE_H
#define CLEAVAGE_SITE_H

#include <stdbool.h>
#include <stddef.h>
#include "anchor-site.h"
#include "fasta.h"
#include "logger.h"
#include "reference.h"

/*
   The start of HA2, which anchors the site, and the bases of the anchor
   taken from the reference.

   The fusion peptide is GLFGAIAG in most haemagglutinins, GIFGAIAG in
   the H3 of horses and pigs, GLLGAIAG in some others and GFFGAIAG in
   influenza B, so the two residues which vary are left free. The
   pattern is found exactly once in each of the 89 HA references of the
   reference sets flucut ships, of every subtype H1 to H16 and of
   influenza B.
*/
#define CLEAVAGE_FUSION_PEPTIDE "G..GAIAG"
#define CLEAVAGE_ANCHOR_NT 24

/* the cleavage site proper: the 8 codons before HA2 */
#define CLEAVAGE_SITE_NT 24
#define CLEAVAGE_SITE_AA (CLEAVAGE_SITE_NT / 3)

/* how a motif reads */
enum { CLEAVAGE_MULTIBASIC,  /* a run of basic residues (R, K) */
       CLEAVAGE_MONOBASIC,   /* one basic residue before HA2 */
       CLEAVAGE_UNREADABLE,  /* a stop codon, an ambiguous base, or no
                                basic residue before HA2 */
       CLEAVAGE_N_READINGS };

/* How the residues <aa> read (CLEAVAGE_MULTIBASIC, ...), and how many
   of their last 6 are basic */
size_t cleavage_reading_of(const char *aa);
size_t cleavage_basic_residues(const char *aa);

/* The name used in the reports */
const char* cleavage_reading_name(size_t reading);

/*
   The site of one HA reference: the anchor which finds it in a read,
   and the site the reference itself carries, which the reads are
   compared with. The reference of a sample is the one closest to its
   reads, so a majority of the reads which does not read like it means
   the sample holds the reads of more than one haemagglutinin.
*/
typedef struct {
  anchor_t anchor;
  char nt[CLEAVAGE_SITE_NT + 1];
  char aa[CLEAVAGE_SITE_AA + 1];
  size_t reading;                  /* CLEAVAGE_MULTIBASIC, ... */
  char fusion[CLEAVAGE_SITE_AA + 1];   /* the anchor, translated */
} cleavage_reference_t;

/* the sites of the HA references of a run, built once and shared */
typedef struct cleavage_anchors_t cleavage_anchors_t;

/*
   Build the anchors of the HA references of <references> (NULL: none),
   which <index> says the segment of; NULL if none of them has the
   fusion peptide (e.g. references which are not whole segments). The
   references must exist as long as the anchors are used.
*/
cleavage_anchors_t* cleavage_anchors_new(const fasta_t *references,
    const reference_index_t *index, Logger *logger);

/*
   The anchor to use for a sample, among the references <allowed> (NULL:
   all of them; see auto_settings_choose_references()): the first HA
   reference allowed which has one, or NULL. With the references of
   several subtypes, one HA reference is allowed per sample, so this is
   the anchor of the subtype of that sample.
*/
const cleavage_reference_t* cleavage_anchors_chosen(
    const cleavage_anchors_t *self, const bool *allowed);

/* How many of the references have an anchor */
size_t cleavage_anchors_n(const cleavage_anchors_t *self);

void cleavage_anchors_delete(cleavage_anchors_t *self);

#endif
