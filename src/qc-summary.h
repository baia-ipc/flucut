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

/* a segment is complete if at least this share of its positions is well
   covered, and low if less than this other share is (e.g. only its ends
   are covered) */
#define QC_COMPLETE_FRACTION 0.9
#define QC_LOW_FRACTION 0.5

/* how well a segment is covered */
typedef enum {
  QC_SEGMENT_MISSING,      /* no reads located in a reference of it */
  QC_SEGMENT_LOW,          /* see QC_LOW_FRACTION */
  QC_SEGMENT_PARTIAL,      /* between low and complete */
  QC_SEGMENT_COMPLETE      /* see QC_COMPLETE_FRACTION */
} qc_segment_status_t;

/* fewer reads than this are too few to tell whether the segments of a
   sample are complete */
#define QC_FEW_READS 1000
/* from this share of the read ends with a primer not located in the
   references, the sample is not one the references describe (e.g. not
   influenza A) */
#define QC_NOT_LOCATED 0.5

/* the state of a sample: that of its segments expected, taken together,
   or the reason why they cannot be judged */
typedef enum {
  QC_SAMPLE_COMPLETE,      /* all of them complete */
  QC_SAMPLE_PARTIAL,       /* all found, some not complete */
  QC_SAMPLE_INCOMPLETE,    /* some missing */
  QC_SAMPLE_FEW_READS,     /* not complete, with fewer than QC_FEW_READS */
  QC_SAMPLE_NOT_LOCATED    /* not complete, see QC_NOT_LOCATED */
} qc_sample_status_t;

/* the size of the name of a subtype (e.g. "H5N1") */
#define QC_SUBTYPE_SIZE 16

/* a segment of a sample */
typedef struct {
  bool found;              /* reads were located in a reference of it */
  qc_segment_status_t status;
  size_t ref;              /* the reference with the greatest depth */
  char subtype[QC_SUBTYPE_SIZE];  /* of that reference ("": not named) */
  size_t ref_len;
  double mean_depth;       /* reads per position */
  size_t median_depth;     /* the depth of half of the positions or more */
  double bases;            /* the bases of the located reads (the sum of
                              the depths) */
  double minus;            /* share of the depth from reads on the minus
                              strand (of those whose strand is known) */
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
  /* the segments expected (all eight, or the only one of a run of one
     segment), those of them which are complete, and the state of the
     sample */
  size_t n_expected, n_complete;
  qc_sample_status_t status;
  /* the share of the read ends with a primer which were not located */
  double not_located;
  int only_segment;        /* 0: all segments are expected */
  /* the subtype of the sample: the H of the reference of its HA and the
     N of that of its NA ("?" for one not known; "": neither) */
  char subtype[QC_SUBTYPE_SIZE];
  /* the references of the segments found name different subtypes (e.g.
     a reassortant, or internal genes closest to another subtype) */
  bool mixed;
} qc_summary_t;

/* The depth of the reference <ref> of <refs> in the counts <c>, and what
   it is made of, as for a segment */
void qc_reference(const fasta_t *refs, const counter_t *c, size_t ref,
    qc_segment_t *out);

/* The summary of the counts <c> for the references <refs> (NULL: none,
   no segment is then found); <only_segment>: the segment of a run of one
   segment alone (1 to 8; 0: all are expected) */
void qc_summary(const fasta_t *refs, const counter_t *c, int only_segment,
    qc_summary_t *out);

/* the shares of the depth of a segment are not given below this mean
   depth: they would be those of a few reads */
#define QC_MIN_SHARE_DEPTH 20

/* Whether the shares of the depth of the segment <q> (internally primed
   reads, internal deletions) mean something */
bool qc_has_shares(const qc_segment_t *q);

/* Whether the segment <segment> (1 to 8) is expected in the summary */
bool qc_expects(const qc_summary_t *self, int segment);

/* The subtype of a sample from the subtypes named by the references of
   its HA (<ha>) and of its NA (<na>), "" if not named, into <buf> (of
   size QC_SUBTYPE_SIZE) */
void qc_subtype(const char *ha, const char *na, char *buf);

/* Unit tests */
void qc_summary_test(void);

#endif
