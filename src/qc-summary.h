/* (c) 2026 Giorgio Gonnella */

/*
   The quality of an influenza A sample, segment by segment, from the
   counts of its reads.

   The trimming locates the reads in the references (amplicon.h): what
   it counts on the way describes the sample itself, whatever is done
   with the trimmed reads afterwards. For each of the eight segments:
   the reference its reads were located in and the subtype named in its
   description; how deep and how completely the segment is covered; and
   how much of that depth comes from reads which do not represent the
   whole segment (internally primed amplicons, genomes with an internal
   deletion). The subtype of the sample is that of the references of its
   HA and NA.

   The coverage is the one reconstructed from the positions of the read
   ends (counter.h), not a base-by-base alignment.
*/

#ifndef QC_SUMMARY_H
#define QC_SUMMARY_H

#include <stdbool.h>
#include <stddef.h>
#include "counter.h"
#include "fasta.h"
#include "segments.h"

/* the depth from which a position is well covered (the reads FLU-minion
   of IRMA needs to support a base of its amended consensus) */
#define QC_MIN_DEPTH 50

/* the size of the name of a subtype (e.g. "H5N1") */
#define QC_SUBTYPE_SIZE 16

/* a segment of a sample */
typedef struct {
  bool found;              /* reads were located in a reference of it */
  size_t ref;              /* the reference with the greatest depth */
  char subtype[QC_SUBTYPE_SIZE];  /* of that reference ("": not named) */
  size_t ref_len;
  double mean_depth;       /* reads per position */
  double covered;          /* share of the positions with a read */
  double well_covered;     /* ... with at least QC_MIN_DEPTH reads */
  double internal;         /* share of the depth from internally primed
                              reads */
  double deletion;         /* ... from reads with an internal deletion */
  size_t complete;         /* complete amplicons (both primers) */
  size_t short_complete;   /* of them, shorter than the reference (see
                              lengths.h) */
} qc_segment_t;

typedef struct {
  qc_segment_t segments[SEGMENTS_N + 1];   /* by segment number (1-8) */
  size_t n_found;          /* segments with located reads */
  /* the subtype of the sample: the H of the reference of its HA and the
     N of that of its NA ("?" for one not known; "": neither) */
  char subtype[QC_SUBTYPE_SIZE];
  /* the references of the segments found name different subtypes (e.g.
     a reassortant, or internal genes closest to another subtype) */
  bool mixed;
} qc_summary_t;

/* The summary of the counts <c> for the references <refs> (NULL: none,
   no segment is then found) */
void qc_summary(const fasta_t *refs, const counter_t *c, qc_summary_t *out);

/* The subtype of a sample from the subtypes named by the references of
   its HA (<ha>) and of its NA (<na>), "" if not named, into <buf> (of
   size QC_SUBTYPE_SIZE) */
void qc_subtype(const char *ha, const char *na, char *buf);

/* Unit tests */
void qc_summary_test(void);

#endif
