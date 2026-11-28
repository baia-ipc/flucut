/* (c) 2026 Giorgio Gonnella */

/*
   What the reads carry at the cleavage site of the haemagglutinin
   (cleavage-site.h), counted.

   What is reported is not an assembly: it is what the reads themselves
   carry at the site, with how many of them agree. It is available
   before any assembly, and tells whether the consensus that follows
   can be trusted there.

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

#ifndef CLEAVAGE_COUNTS_H
#define CLEAVAGE_COUNTS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include "cleavage-site.h"

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

/* the motifs listed in the reports, most frequent first */
#define CLEAVAGE_TOP_MOTIFS 5

/*
   What the reads must reach for a reading of the site (the residues) to
   count as a site the sample holds rather than an error of a read: the
   share of the reads that carry it, and at least this many reads. Every
   misread read has a reading of its own, so the number of readings seen
   is mostly the number of errors; the readings above this share are the
   ones worth looking at.

   Chosen on 104 Nanopore samples (docs/cleavage-site.md). Errors of one
   residue recur in the same samples at 5 to 13% of the reads, so a
   share of 1% or 5% still counts them (86 and 20 of the 104 samples
   have more than 4 and more than 2 readings), while at 10% most samples
   have one reading and the others are mostly those holding a second
   haemagglutinin.
*/
#define CLEAVAGE_MIN_SITE_SHARE 0.10
#define CLEAVAGE_MIN_SITE_READS 2

/* how far the motif reported can be trusted */
enum { CLEAVAGE_FIRM,       /* enough reads agree */
       CLEAVAGE_CONFIRM,    /* the motif is the majority, but too few of
                               the reads agree for it to stand alone */
       CLEAVAGE_TOO_FEW_READS,
       CLEAVAGE_N_VERDICTS };

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
  /* distinct spellings, and distinct residues, among them: both include
     the readings of single misread reads, so they say how noisy the
     reads are, not how many haemagglutinins they hold */
  size_t distinct;
  size_t distinct_proteins;
  /* the residues carried by at least CLEAVAGE_MIN_SITE_SHARE of the
     reads and CLEAVAGE_MIN_SITE_READS reads: the sites the sample holds */
  size_t supported_proteins;
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
const char* cleavage_verdict_name(size_t verdict);

#endif
