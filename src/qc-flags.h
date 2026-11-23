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
#include "counter.h"
#include "qc-summary.h"

/* at least this share of the reads without a primer (e.g. the scheme is
   not that of the run, or the reads were trimmed already) */
#define QC_FLAGS_NO_PRIMER 0.2
/* at least this share of the depth of a segment from internally primed
   reads, from reads with an internal deletion */
#define QC_FLAGS_INTERNAL 0.1
#define QC_FLAGS_DELETION 0.5

#define QC_FLAGS_MAX 8
#define QC_FLAGS_TEXT_SIZE 160

typedef struct {
  const char *key;                 /* a name for programs, e.g. few_reads */
  char text[QC_FLAGS_TEXT_SIZE];   /* for people */
} qc_flag_t;

typedef struct {
  qc_flag_t flags[QC_FLAGS_MAX];
  size_t n;
} qc_flags_t;

/* The flags of a sample with the counts <s> and the summary <qc>
   (<refs>: with references, otherwise only the flags of the counts) */
void qc_flags(const counter_reads_t *s, const qc_summary_t *qc, bool refs,
    qc_flags_t *out);

/* Unit tests */
void qc_flags_test(void);

#endif
