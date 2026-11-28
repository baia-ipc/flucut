/* (c) 2026 Giorgio Gonnella */

/*
   What the reads carry at each marker (markers.h), counted.

   For every marker the reads are counted by the value they give: the
   residue at the position, the residues at the end of a protein, the
   length of a reading frame or of a stretch. A value is kept with the
   reads which gave it, so the result is a distribution and not a call:
   a residue carried by a part of the reads is reported with its share,
   which is where this differs from reading the same position off a
   consensus, where such a part disappears.

   Beside each value the reads are counted by the residues they carry
   where the anchor matched. The anchor is conserved, so the share of
   the reads agreeing there is how much the reads of this sample agree
   where nothing is difficult: the floor below which a value carried by
   few reads is not told from the errors of the reads.

   The memory does not depend on the number of reads: once
   MARKER_VALUES values are known, a read with another one is counted
   as a read but not as a new value (<overflow>), which can only lower
   the shares reported, never raise them.
*/

#ifndef MARKER_COUNTS_H
#define MARKER_COUNTS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include "markers.h"

/* the distinct values kept, and the slots of the table holding them */
#define MARKER_VALUES 128
#define MARKER_SLOTS 256

/* what one window counted */
typedef struct {
  size_t reads;
  size_t overflow;
  size_t n;
  struct {
    char value[MARKER_VALUE_MAX + 1];
    size_t reads;
  } entries[MARKER_VALUES];
  size_t slots[MARKER_SLOTS];   /* entry + 1, 0 empty */
} marker_tally_t;

/* what the reads gave for one marker */
typedef struct {
  size_t searched;            /* reads the marker was looked for in */
  marker_tally_t value;
  marker_tally_t anchor;      /* the anchor as the reads carry it */
  /* what the reference of the sample carries ("": not known) */
  char reference_value[MARKER_VALUE_MAX + 1];
  char reference_anchor[MARKER_ANCHOR_AA + 1];
} marker_one_counts_t;

/*
   The markers read from the reads of a group of input files. A counts
   belongs to one thread; the counts of the threads are summed
   (marker_counts_merge()) in the order of the input.
*/
typedef struct {
  marker_one_counts_t markers[MARKER_N_MARKERS];
} marker_counts_t;

void marker_counts_clear(marker_counts_t *self);

/* Record what the reference the reads are read against carries for
   <marker>; ignored if it is already recorded */
void marker_counts_set_reference(marker_counts_t *self, size_t marker,
    const char *value, const char *anchor_aa);

/* Count a read <marker> was looked for in, whether or not it gave a
   value */
void marker_counts_searched(marker_counts_t *self, size_t marker);

/* Count the value <value> of <marker> read from a read, and the
   residues <aa> the read carries where its anchor matched */
void marker_counts_add(marker_counts_t *self, size_t marker,
    const char *value);
void marker_counts_add_anchor(marker_counts_t *self, size_t marker,
    const char *aa);

/* The reads <marker> was read from (the values and the overflow) */
size_t marker_counts_reads(const marker_counts_t *self, size_t marker);

/* Whether any read was searched for any marker */
bool marker_counts_any(const marker_counts_t *self);

void marker_counts_merge(marker_counts_t *self, const marker_counts_t *other);

/* Save to <file>, and load one of its lines (as counter_save(), for the
   resumption of a run in progress); false if the line is not one */
void marker_counts_save(const marker_counts_t *self, FILE *file);
bool marker_counts_load_line(marker_counts_t *self, char *line);

/* --- what the counts say --- */

/* fewer reads than this: no value is reported */
#define MARKER_MIN_READS 20

/*
   What a value must reach to count as one the sample holds rather than
   an error of a read: the share of the reads carrying it, and at least
   this many reads. Every misread read has a value of its own, so most
   of the values seen are errors; those above this share are the ones
   worth looking at. A minority below it is still reported in the
   distribution, as a value, but not as one the sample holds.

   A tenth, as at the cleavage site: below it a value is not told from
   the disagreement the reads of a sample already show at the anchor,
   where nothing is difficult. On 104 Nanopore samples that
   disagreement is a median of 14% of the reads on R10 and 39% on R9
   (docs/markers.md).
*/
#define MARKER_MIN_SHARE 0.10
#define MARKER_MIN_VALUE_READS 3

/*
   What the reads must reach for the value reported to be firm: they
   must agree on it, and they must agree where the anchor matched,
   which is where nothing is difficult. A sample failing the first has
   a position its reads cannot resolve; one failing the second has
   reads too noisy, or too mixed, for a claim about any position.

   Below MARKER_MIN_MAJORITY the most frequent value is not a majority
   of anything: no value is reported. A window of many codons reaches
   that state on noisy reads, because one indel anywhere in it changes
   what is read, so the markers read over a long window are often
   unresolved on Nanopore reads and resolved on Illumina ones
   (docs/markers.md).
*/
#define MARKER_MIN_SUPPORT 0.60
#define MARKER_MIN_ANCHOR_SUPPORT 0.60
#define MARKER_MIN_MAJORITY 0.35

/* how far the value reported can be trusted */
enum { MARKER_FIRM,       /* enough reads agree */
       MARKER_CONFIRM,    /* the value is the majority, but too few of
                             the reads agree for it to stand alone */
       MARKER_UNRESOLVED, /* the reads do not agree on any value: the
                             one reported is not the sample's */
       MARKER_TOO_FEW_READS,
       MARKER_N_VERDICTS };

/* the values listed in the reports, most frequent first */
#define MARKER_TOP_VALUES 4

typedef struct {
  char value[MARKER_VALUE_MAX + 1];
  size_t reads;
  double share;
  bool attention;      /* the value deserves attention (markers.h) */
} marker_value_t;

/* what the reads say about one marker */
typedef struct {
  bool found;                    /* false: no read gave a value */
  size_t searched;
  size_t reads;                  /* reads it was read from */
  size_t distinct;               /* distinct values among them, the ones
                                    of single misread reads included */
  size_t supported;              /* values above MARKER_MIN_SHARE: the
                                    values the sample holds, more than
                                    one of them a mixture */
  char value[MARKER_VALUE_MAX + 1];   /* the most frequent */
  size_t value_reads;
  double support;
  marker_value_t top[MARKER_TOP_VALUES];
  size_t n_top;
  /* the residues the reads carry where the anchor matched, the most
     frequent of them, and the share of the reads agreeing: the noise
     floor of this sample at this marker. The residues of the anchor in
     the reference are beside them, because an anchor which matched
     somewhere else gives reads which agree with each other and not
     with the reference */
  size_t anchor_reads;
  double anchor_support;
  bool anchor_weak;
  char anchor_value[MARKER_ANCHOR_AA + 1];
  char reference_anchor[MARKER_ANCHOR_AA + 1];
  bool anchor_differs;
  /* what the reference carries, and whether the most frequent value
     differs from it */
  char reference_value[MARKER_VALUE_MAX + 1];
  bool differs_from_reference;
  /* a value the sample holds deserves attention, and the share of the
     reads carrying such a value */
  bool attention;
  double attention_share;
  size_t verdict;                /* MARKER_FIRM, ... */
} marker_result_t;

/* What <counts> says about <marker> */
void marker_result_of(const marker_counts_t *counts, size_t marker,
    marker_result_t *out);

/* The names used in the reports */
const char* marker_verdict_name(size_t verdict);

#endif
