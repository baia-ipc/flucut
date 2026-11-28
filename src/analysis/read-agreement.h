/* (c) 2026 Giorgio Gonnella */

/*
   The agreement of the probes of a read, in the choice of the reference
   of each segment (see auto_settings_choose_references): how often the
   probes of one read vote for the same reference.

   A segment whose votes are split between two references has two
   possible causes, which the votes alone do not tell apart:
   - one strain lying between the two references: the probes of the same
     read go to different references;
   - two strains in the sample: the probes of a read agree, the reads
     disagree with each other.
   For each segment, the reads with several voting probes are counted,
   and those whose probes all vote for the same reference, by reference.
   These are compared to the numbers expected if each probe voted
   independently of its read, with the vote shares of the segment (a
   read with k voting probes is then for a reference of share p with
   probability p^k). With one strain the observed numbers are about as
   expected; with two strains they are higher, and the excess of the
   reads which agree on the second reference, as a share of the reads,
   estimates the fraction of the minor strain.

   The statistic is reported, not used: it does not change the
   references chosen.
*/

#ifndef READ_AGREEMENT_H
#define READ_AGREEMENT_H

#include <stddef.h>
#include "fasta.h"
#include "logger.h"
#include "segments.h"

/* the most voting probes of a read which are considered */
#define READ_AGREEMENT_MAX_PROBES 16

typedef struct {
  size_t n_refs;
  /* for each reference, the reads with several voting probes in its
     segment which all vote for it (owned) */
  size_t *agreeing;
  /* for each segment, the reads by number of voting probes in it (from
     2) */
  size_t reads[SEGMENTS_N + 1][READ_AGREEMENT_MAX_PROBES + 1];
} read_agreement_t;

void read_agreement_init(read_agreement_t *self, size_t n_refs);

/* Add a read: <read_votes> are the references its <n> voting probes
   voted for (those beyond READ_AGREEMENT_MAX_PROBES are not considered),
   <segment> the segment of each reference */
void read_agreement_add(read_agreement_t *self, const int *segment,
    const size_t *read_votes, size_t n);

/* The reads with several voting probes in the segment <s> */
size_t read_agreement_reads(const read_agreement_t *self, int s);

/* The reads of the segment <s> expected to agree on a reference which
   has the share <share> of the votes of the segment, if the probes
   voted independently of their read */
double read_agreement_expected(const read_agreement_t *self, int s,
    double share);

/* Report the agreement in each segment with several references, for the
   reads of <file>: <segment> and <votes> have a value for each of the
   references <fa> */
void read_agreement_log(const read_agreement_t *self, const fasta_t *fa,
    const int *segment, const size_t *votes, const char *file,
    Logger *logger);

void read_agreement_free(read_agreement_t *self);

/* Unit tests */
void read_agreement_test(void);

#endif
