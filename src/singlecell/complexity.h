/* (c) 2026 Giorgio Gonnella */

/*
   Library complexity of a sample: estimation of the number of distinct
   molecules of the library, from the number of reads and of distinct
   molecule keys (e.g. the UMI and the first and last bases of the
   insert, which identify the copies of a molecule without an alignment),
   for deciding whether deeper sequencing of a library is worthwhile.

   Sequencing errors make copies of a molecule look distinct: a key of
   length L is without errors with probability p = (1 - e)^L (e: error
   rate per base, e.g. estimated from the adaptor alignments); each read
   with an error in the key is assumed to add a distinct key. Thus, the
   pN reads with correct keys contain D = D_obs - (1 - p)N distinct
   molecules. The number of molecules M is then estimated as in Picard
   (Lander-Waterman equation, sampling with replacement):
   D = M (1 - exp(-pN / M)).

   The estimate is approximate: it assumes molecules sampled uniformly
   (no amplification bias) and independent errors; with few duplicates,
   it is sensitive to the error of the count of the distinct keys. With
   few keys without errors (p < COMPLEXITY_MIN_ERROR_FREE, e.g. more than
   about 5% errors with keys of 24 bases), the correction dominates the
   count and there is no estimate.
*/

#ifndef COMPLEXITY_H
#define COMPLEXITY_H

#include <stdbool.h>
#include <stddef.h>

/* minimum number of reads with correct keys for an estimate */
#define COMPLEXITY_MIN_READS 100

/* minimum fraction of the keys without errors for an estimate */
#define COMPLEXITY_MIN_ERROR_FREE 0.3

typedef struct {
  bool estimated;          /* false: too few reads or too many errors */
  bool high_error_rate;    /* not estimated: p < COMPLEXITY_MIN_ERROR_FREE
                              (with at least COMPLEXITY_MIN_READS reads) */
  double error_free;       /* p: fraction of the keys without errors */
  double molecules_seen;   /* distinct molecules in the N reads */
  double duplication;      /* fraction of duplicate reads */
  bool saturated_unknown;  /* no duplicates: M cannot be estimated */
  double molecules;        /* M: molecules of the library */
  double molecules_2x;     /* distinct molecules expected with 2N reads */
} complexity_t;

/*
   Estimate the complexity from <reads> reads with <distinct_keys>
   distinct keys, of mean length <mean_key_length>, and the error rate per
   base <error_rate>
*/
void complexity_estimate(size_t reads, double distinct_keys,
    double mean_key_length, double error_rate, complexity_t *out);

/* Unit tests */
void complexity_test(void);

#endif
