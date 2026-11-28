/* (c) 2026 Giorgio Gonnella */

/*
   What deserves attention in a sample: short remarks derived from its
   counts and from its segments (qc-summary.h), for the reports. A flag
   does not change what is written: it points the reader to the numbers
   worth a look (e.g. a segment whose depth is mostly internally primed
   reads). The thresholds are those below and those of the state of a
   sample (QC_FEW_READS, QC_NOT_LOCATED).
*/

#ifndef QC_FLAGS_H
#define QC_FLAGS_H

#include <stddef.h>
#include "cleavage-counts.h"
#include "counter.h"
#include "qc-summary.h"
#include "subtype-reads.h"

/* at least this share of the reads without a primer (e.g. the scheme is
   not that of the run, or the reads were trimmed already) */
#define QC_FLAGS_NO_PRIMER 0.2
/* at least this share of the depth of a segment from internally primed
   reads, from reads with an internal deletion */
#define QC_FLAGS_INTERNAL 0.1
#define QC_FLAGS_DELETION 0.5
/* below this share of the reads placed on a segment of the references
   (counter.h), the references describe only a small part of the sample.
   Chosen from the 104 Nanopore samples of the multisegment RT-PCR of the
   HA I351K validation: the share has median 0.992 and is 0.63 or more in
   102 of them; it is 0.014 and 0.154 in the two samples whose reads map
   to influenza at 1.3% and 15% (the other 102: 63% or more). The
   threshold lies at about twice the second and half the third */
#define QC_FLAGS_PLACED 0.3

#define QC_FLAGS_MAX 8
#define QC_FLAGS_TEXT_SIZE 160

typedef struct {
  const char *key;                 /* a name for programs, e.g. few_reads */
  char text[QC_FLAGS_TEXT_SIZE];   /* for people */
  /* the results of the sample rest on a small part of its reads: shown
     as a problem, not only as a remark */
  bool severe;
} qc_flag_t;

typedef struct {
  qc_flag_t flags[QC_FLAGS_MAX];
  size_t n;
} qc_flags_t;

/*
   The flags of a sample with the counts <s> and the summary <qc>
   (<refs>: with references, otherwise only the flags of the counts).

   <site> is what the reads said about the cleavage site of the
   haemagglutinin (cleavage-counts.h; NULL: not read): a site which is not
   firm is pointed out here as well, so that it is seen in the run
   overview and not only in the sample.

   <subtypes> are the reads of each HA and NA subtype among the first
   reads of the sample (subtype-reads.h; NULL: not counted): reads of a
   second subtype above the threshold are pointed out, as what was
   observed, since their cause is not known.
*/
void qc_flags(const counter_reads_t *s, const qc_summary_t *qc, bool refs,
    const cleavage_result_t *site, const subtype_reads_t *subtypes,
    qc_flags_t *out);

/* The share of the reads analyzed (after splitting) placed on a segment
   of the references; 0 without reads */
double qc_flags_reads_placed(const counter_reads_t *s);

/* Unit tests */
void qc_flags_test(void);

#endif
