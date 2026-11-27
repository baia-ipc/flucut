/* (c) 2026 Giorgio Gonnella */

/*
   The cleavage site of the haemagglutinin, read from the reads.

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

   What is reported here is not an assembly: it is what the reads
   themselves carry at the site, with how many of them agree. It is
   available before any assembly, and tells whether the consensus that
   follows can be trusted there.

   Where it is read from
   ---------------------

   HA2 begins directly after the cleavage site with the fusion peptide
   (GLFGAIAG and a few variants, see CLEAVAGE_FUSION_PEPTIDE), which
   every haemagglutinin has and which lies outside the low-complexity
   stretch. Those bases of the HA
   reference are the anchor (anchor-site.h); the cleavage site is the
   CLEAVAGE_SITE_NT bases of the read immediately before it, read back
   in whole codons from the first base of HA2, so the frame is set by
   the anchor and no indel elsewhere in the read can shift it.

   Which reads are read is left to the caller, because it is not the
   same on the two platforms (see trimmer.c): the reads of whole
   amplicons (Nanopore) reach the site only when they are full
   templates, as the internally primed ones begin just after the site
   and the genomes with an internal deletion have removed it, while a
   read of a fragment (Illumina) which covers the site usually carries
   no primer at all and has no category, so there every read is
   searched and the anchor itself says which ones cover the site.

   How far it can be trusted
   -------------------------

   Two numbers decide: how many reads gave the site, and what share of
   them agree on the same motif (its support). How much support is
   enough depends on how often the reads make an indel there, which
   differs between Illumina, R10 and R9 reads, so each is judged by its
   own number (see the thresholds below and docs/cleavage-site.md).
*/

#ifndef CLEAVAGE_SITE_H
#define CLEAVAGE_SITE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
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

/*
   The control windows: two of the width of the site, taken from the
   same reads on the other side of the anchor, where HA2 is ordinary
   coding sequence. They say how much the reads of this sample agree
   where nothing is difficult, which tells a site the reads cannot
   resolve from reads that agree nowhere.

   They are reported, not used to judge the site. On 104 samples the
   site agreed less than the control in nearly every one, in the
   samples whose consensus was right as much as in those whose consensus
   was wrong, so the excess does not separate them and the ratio of the
   two is a worse test than the agreement at the site itself
   (docs/cleavage-site.md).
*/
#define CLEAVAGE_N_CONTROLS 2
#define CLEAVAGE_CONTROL_MOTIFS 256
#define CLEAVAGE_CONTROL_SLOTS 512

/* what one window of residues counted: its reads by the residues they
   carry, bounded as the table of the site is. Used for the control
   windows and for the fusion peptide */
typedef struct {
  size_t reads;
  size_t overflow;
  size_t n;
  struct {
    char aa[CLEAVAGE_SITE_AA + 1];
    size_t reads;
  } entries[CLEAVAGE_CONTROL_MOTIFS];
  size_t slots[CLEAVAGE_CONTROL_SLOTS];
} cleavage_tally_t;

/*
   The distinct motifs kept, and the slots of the table that holds them.
   The memory does not depend on the number of reads: once
   CLEAVAGE_MAX_MOTIFS motifs are known, a read with another one is
   counted as a read but not as a new motif (cleavage_counts_t.overflow),
   which can only lower the support reported, never raise it. The motif
   of the sample is far more frequent than one in CLEAVAGE_MAX_MOTIFS,
   so it is always among those kept.
*/
#define CLEAVAGE_MAX_MOTIFS 1024
#define CLEAVAGE_SLOTS 2048

/* fewer reads than this: no motif is reported */
#define CLEAVAGE_MIN_READS 20

/*
   What a result must reach before it is reported as firm: the reads
   must agree at the site, and they must agree where nothing is
   difficult. The second is what the control windows are for.

   The two are independent. A sample whose reads agree in the control
   windows but not at the site has a site its reads cannot resolve; a
   sample whose reads agree in neither has reads too noisy, or too
   mixed, for any claim about any position, and the site is then not
   wrong so much as unknown.

   Neither number depends on the platform or the chemistry. On 104
   Nanopore samples the control windows separate the chemistries by
   themselves - every R9 sample falls below the floor and the R10
   samples sit well above it - so the data says what a rule naming the
   chemistry used to say, and says it from the reads in hand
   (docs/cleavage-site.md).
*/
#define CLEAVAGE_MIN_SUPPORT 0.45
#define CLEAVAGE_MIN_CONTROL_SUPPORT 0.60

/* the kind of data, which decides the support needed */
enum { CLEAVAGE_CHEMISTRY_UNKNOWN,  /* Nanopore, chemistry not in the titles */
       CLEAVAGE_CHEMISTRY_R9,
       CLEAVAGE_CHEMISTRY_R10,
       CLEAVAGE_CHEMISTRY_ILLUMINA,
       CLEAVAGE_N_CHEMISTRIES };

/* how a motif reads */
enum { CLEAVAGE_MULTIBASIC,  /* a run of basic residues (R, K) */
       CLEAVAGE_MONOBASIC,   /* one basic residue before HA2 */
       CLEAVAGE_UNREADABLE,  /* a stop codon, an ambiguous base, or no
                                basic residue before HA2 */
       CLEAVAGE_N_READINGS };

/* the motifs listed in the reports, most frequent first */
#define CLEAVAGE_TOP_MOTIFS 5

/* how far the motif reported can be trusted */
enum { CLEAVAGE_FIRM,       /* enough reads agree */
       CLEAVAGE_CONFIRM,    /* the motif is the majority, but too few of
                               the reads agree for it to stand alone */
       CLEAVAGE_TOO_FEW_READS,
       CLEAVAGE_N_VERDICTS };

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

/*
   The motifs read from the reads of a group of input files. A counts
   belongs to one thread; the counts of the threads are summed
   (cleavage_counts_merge()) in the order of the input.
*/
typedef struct {
  size_t searched;   /* reads the site was looked for in */
  size_t overflow;   /* reads whose motif found no slot (see above) */
  /* the control windows, each counted by residues in a table of its
     own; smaller than the one of the site, which is enough to hold a
     majority, and bounded in the same way */
  cleavage_tally_t control[CLEAVAGE_N_CONTROLS];
  /* the fusion peptide itself, as the reads carry it: the anchor is
     conserved, so this is a check that it landed where it should and
     shows the variant this sample carries */
  cleavage_tally_t fusion;
  /* the fusion peptide of the HA reference ("": none) */
  char reference_fusion[CLEAVAGE_SITE_AA + 1];
  size_t n_motifs;
  /* the site of the HA reference the reads were read against ("": none) */
  char reference_aa[CLEAVAGE_SITE_AA + 1];
  struct {
    char nt[CLEAVAGE_SITE_NT + 1];
    size_t reads;
  } motifs[CLEAVAGE_MAX_MOTIFS];
  size_t slots[CLEAVAGE_SLOTS];   /* motif + 1, 0 empty */
} cleavage_counts_t;

void cleavage_counts_clear(cleavage_counts_t *self);

/*
   Record the site of the HA reference the reads are being read against
   (<aa>) and its fusion peptide (<fusion>, the anchor); ignored if one
   is already recorded
*/
void cleavage_counts_set_reference(cleavage_counts_t *self, const char *aa,
    const char *fusion);

/* Count a read searched, whether or not the site was found in it */
void cleavage_counts_searched(cleavage_counts_t *self);

/* Count the motif <nt> (CLEAVAGE_SITE_NT bases) read from a read */
void cleavage_counts_add(cleavage_counts_t *self, const char *nt);

/* Count the bases <nt> read from the control window <i> of a read */
void cleavage_counts_add_control(cleavage_counts_t *self, size_t i,
    const char *nt);

/* Count the bases <nt> the read carries where the anchor matched */
void cleavage_counts_add_fusion(cleavage_counts_t *self, const char *nt);

/* The reads the site was read from (the motifs and the overflow) */
size_t cleavage_counts_reads(const cleavage_counts_t *self);

void cleavage_counts_merge(cleavage_counts_t *self,
    const cleavage_counts_t *other);

/* Save to <file>, and load one of its lines (as counter_save(), for the
   resumption of a run in progress); false if the line is not one */
void cleavage_counts_save(const cleavage_counts_t *self, FILE *file);
bool cleavage_counts_load_line(cleavage_counts_t *self, char *line);

/* one of the sites the reads carry */
typedef struct {
  char aa[CLEAVAGE_SITE_AA + 1];   /* the residues */
  char nt[CLEAVAGE_SITE_NT + 1];   /* their most frequent spelling */
  size_t reads;                    /* reads whose site reads as this */
  size_t reading;                  /* CLEAVAGE_MULTIBASIC, ... */
} cleavage_motif_t;

/*
   What the reads say about the site.

   The reads are counted by the residues they carry, not by the bases:
   the residues are what decides the pathotype, and the reads agree on
   them far more often than on the spelling, which differs by a base
   wherever a read has an error. <nt> is the most frequent spelling of
   the residues reported.
*/
typedef struct {
  bool found;                      /* false: no read gave the site */
  size_t searched;
  size_t reads;                    /* reads the site was read from */
  size_t distinct;                 /* distinct spellings among them */
  size_t distinct_proteins;        /* distinct residues among them */
  char aa[CLEAVAGE_SITE_AA + 1];   /* the most frequent residues */
  char nt[CLEAVAGE_SITE_NT + 1];   /* their most frequent spelling */
  size_t motif_reads;              /* reads carrying those residues */
  double support;                  /* motif_reads / reads */
  double multibasic_share;         /* reads whose own site is multibasic */
  /* the control windows, pooled: the reads they were read from and the
     share of those agreeing on their most frequent residues. Reported
     beside the site, as what the reads of this sample agree on where
     nothing is difficult; not used to judge the site */
  size_t control_reads;
  double control_support;          /* 0 where there are no control reads */
  /* the reads of this sample do not agree even in the control windows:
     too noisy, or holding more than one virus, for a claim about any
     position */
  bool control_weak;
  size_t basic;                    /* basic residues in the last 6 of it */
  size_t reading;                  /* CLEAVAGE_* */
  size_t verdict;                  /* CLEAVAGE_FIRM, ... */
  size_t chemistry;                /* the chemistry it was judged by */
  /* the site of the HA reference of the sample ("": not known), and
     whether the majority of the reads reads like it (multibasic against
     multibasic, monobasic against monobasic). It does not when the
     sample holds the reads of more than one haemagglutinin, and the
     motif reported is then not the one of its own virus, so such a
     result is never firm. */
  char reference_aa[CLEAVAGE_SITE_AA + 1];
  size_t reference_reading;
  bool disagrees_with_reference;
  /* the fusion peptide: the one of the reference, which is the anchor,
     and the one the reads carry where it matched, with the share of
     them agreeing. The two differ when the sample is not the subtype
     of its reference; a low share means the anchor is matching reads
     it should not */
  char reference_fusion[CLEAVAGE_SITE_AA + 1];
  char fusion_aa[CLEAVAGE_SITE_AA + 1];
  size_t fusion_reads;
  double fusion_support;
  /* the most frequent sites, most frequent first: a sample with the
     reads of more than one haemagglutinin in it (a co-infection, or
     reads of another barcode) shows here, as the fusion peptide which
     anchors the site is the same in every subtype */
  cleavage_motif_t top[CLEAVAGE_TOP_MOTIFS];
  size_t n_top;
} cleavage_result_t;

/*
   What <counts> says, judged as data of <chemistry>
   (CLEAVAGE_CHEMISTRY_*)
*/
void cleavage_result_of(const cleavage_counts_t *counts, size_t chemistry,
    cleavage_result_t *out);

/*
   The kind of data: Illumina if <illumina>, otherwise the chemistry
   named by the basecalling model found in the titles of the reads
   (<basecall_model>, NULL or "" if there is none, see read-metadata.h):
   a model of the form dna_r10.4.1_... or dna_r9.4.1_...
*/
size_t cleavage_chemistry(bool illumina, const char *basecall_model);

/* The names used in the reports */
const char* cleavage_chemistry_name(size_t chemistry);
const char* cleavage_reading_name(size_t reading);
const char* cleavage_verdict_name(size_t verdict);

/* Unit tests */
void cleavage_site_test(void);

#endif
