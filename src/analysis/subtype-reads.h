/* (c) 2026 Giorgio Gonnella */

/*
   The reads of a second HA or NA subtype in a sample of influenza A.

   The reference of each segment is chosen by the votes of probes taken
   from the first reads of a sample (auto-settings.h). Here the probes
   of each read which vote in the HA (or in the NA) are looked at
   together: a read with several of them is a read of a subtype if all
   of them vote for references of that subtype (H1 to H19 for the HA, N1
   to N11 for the NA; the references of one subtype are counted
   together, whatever the other half of their name: the HA of an H5N1
   and of an H5N8 reference are both H5).

   The subtype with most votes is the leading one. Reads of another
   subtype have two causes:
   - errors and conserved stretches: a probe now and then matches a
     reference of another subtype best. If each probe did so
     independently of its read, a read with k voting probes would be of
     a subtype with the share p of the votes with probability p^k,
     which gives the reads expected;
   - reads of a second virus in the sample: their probes agree, and the
     reads observed exceed the expected ones.
   The excess (observed minus expected) of the subtype where it is
   largest, the second subtype, is the evidence reported. It is above
   the threshold when it is at least SUBTYPE_READS_MIN_EXCESS reads and
   at least SUBTYPE_READS_MIN_SHARE of the reads with several voting
   probes (see below).

   What it does not say: where the reads of the second subtype come
   from. A mixed infection, a contamination of the sample and barcode
   crosstalk between the samples of a run all give the same reads. Two
   viruses of the same subtype are not told apart, and the other six
   segments are not looked at: their references differ too little
   between subtypes for the probes of a read to agree on one.

   A read whose probes do not agree on any subtype is counted apart.
   Such reads are many when no reference of the set is close to the
   segment: the probes of one read then match here one subtype and
   there another, and the subtype named for that segment rests on few
   reads. That is reported as a distant reference
   (SUBTYPE_READS_SCATTERED_SHARE below), and it is what separates a
   divergent segment from a second virus: a second virus gives reads
   which agree among themselves on another subtype, a distant reference
   gives reads which agree on nothing.

   The counts are those of a sample (subtype_reads_t), kept with its
   other counts (counter.h); the evidence is computed from them when
   the reports are written. Nothing here changes the references chosen.
*/

#ifndef SUBTYPE_READS_H
#define SUBTYPE_READS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include "fasta.h"
#include "logger.h"

/* the two segments whose reference names a subtype of its own */
enum { SUBTYPE_READS_HA, SUBTYPE_READS_NA, SUBTYPE_READS_N_SEGMENTS };

/* the highest number of a subtype (H19; the NA has fewer) */
#define SUBTYPE_READS_MAX_NUMBER 19
/* the most voting probes of a read which are considered */
#define SUBTYPE_READS_MAX_PROBES 16

/*
   The thresholds of the excess of the second subtype
   (benchmark/RESULTS.md, "Reads of a second HA or NA subtype").

   The fewest reads: 3 where none is expected. Simulated samples of one
   subtype had at most 1 read of another one, and in sequenced runs
   which held samples of several subtypes 1 or 2 such reads were common
   (0.4 to 0.7% of the reads counted).

   Their smallest share of the reads with several voting probes: it
   keeps the rule the same when more reads are counted (up to 20,000
   instead of 2000, see auto-settings.h), where that background alone
   gives more than 3 reads.
*/
#define SUBTYPE_READS_MIN_EXCESS 2.5
#define SUBTYPE_READS_MIN_SHARE 0.01

/*
   When the reference of a segment is called distant: the share of the
   reads with several voting probes whose probes agree on no subtype,
   and the fewest such reads for the share to mean anything.

   Measured on 177 HA and NA segments of 104 sequenced samples, each
   compared with the identity of its IRMA consensus to the closest
   record of the set (subtype_check/ of the analyses): no segment with
   a record at least 90% identical reached 0.15, and every segment
   above 0.5 had its closest record at most 88% identical. The
   converse does not hold: most segments with a distant record stay
   below 0.5, so this says that the reads of a segment fit no record
   of the set, not how distant the closest one is.
*/
#define SUBTYPE_READS_SCATTERED_SHARE 0.5
#define SUBTYPE_READS_SCATTERED_MIN_READS 20

/* the counts of the HA or of the NA; the subtypes by their number */
typedef struct {
  /* the voting probes located in the references of each subtype */
  size_t votes[SUBTYPE_READS_MAX_NUMBER + 1];
  /* the reads with several voting probes, all for that subtype */
  size_t agreeing[SUBTYPE_READS_MAX_NUMBER + 1];
  /* the reads by number of voting probes (from 2) */
  size_t reads[SUBTYPE_READS_MAX_PROBES + 1];
} subtype_reads_segment_t;

typedef struct {
  subtype_reads_segment_t segments[SUBTYPE_READS_N_SEGMENTS];
  /* the samples counted: the evidence is that of one sample (the sum
     of several, as in the counts of a whole run, says nothing: their
     subtypes may differ) */
  size_t samples;
} subtype_reads_t;

/* what the counts of the HA or of the NA of a sample say */
typedef struct {
  /* false: nothing to say (no read with several voting probes, or the
     counts are not those of one sample) */
  bool known;
  int leading, second;     /* numbers of the subtypes; 0: no second one */
  size_t reads;            /* the reads with several voting probes */
  size_t leading_reads, second_reads;   /* ... all for that subtype */
  /* the reads whose probes agree on no subtype, and their share of
     <reads>: many of them mean that no reference is close */
  size_t scattered_reads;
  double scattered_share;
  bool distant_reference;
  /* the reads of the second subtype expected if the probes voted
     independently of their read, and the observed ones beyond them */
  double second_expected, second_excess;
  double second_share;     /* the excess as a share of <reads> */
  bool above_threshold;
} subtype_reads_evidence_t;

void subtype_reads_clear(subtype_reads_t *self);

/* The number of the subtype which a reference of the segment <segment>
   with the FASTA description <description> stands for: the H of an HA,
   the N of an NA of influenza A; 0 for any other reference */
int subtype_reads_number(const char *description, int segment);

/* The numbers of all the references <fa>, whose segments are <segment>
   (a new array, to be freed by the caller) */
int* subtype_reads_numbers(const fasta_t *fa, const int *segment);

/* Add a read: <read_votes> are the references its <n> voting probes
   voted for, <segment> and <number> the segment and the subtype number
   of each reference */
void subtype_reads_add(subtype_reads_t *self, const int *segment,
    const int *number, const size_t *read_votes, size_t n);

/* Add the counts of <other> */
void subtype_reads_merge(subtype_reads_t *self,
    const subtype_reads_t *other);

/* The evidence in the segment <which> (SUBTYPE_READS_HA or _NA) */
void subtype_reads_evidence(const subtype_reads_t *self, size_t which,
    subtype_reads_evidence_t *out);

/* Whether the excess <excess> among <reads> reads with several voting
   probes is above the threshold */
bool subtype_reads_above_threshold(double excess, size_t reads);

/* Whether <scattered> of <reads> reads whose probes agree on no
   subtype are enough to call the reference of the segment distant */
bool subtype_reads_distant_reference(size_t scattered, size_t reads);

/* The gene (HA, NA) and the letter of the subtypes (H, N) of <which> */
const char* subtype_reads_gene(size_t which);
char subtype_reads_letter(size_t which);

/* Report the evidence in the reads of <file> */
void subtype_reads_log(const subtype_reads_t *self, const char *file,
    Logger *logger);

/* Save the counts to <file> (text lines starting with "subtype"), and
   add those of such a line to <self>; false if the line is not valid */
void subtype_reads_save(const subtype_reads_t *self, FILE *file);
bool subtype_reads_load_line(subtype_reads_t *self, const char *line);

/* Unit tests */
void subtype_reads_test(void);

#endif
