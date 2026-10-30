/* (c) 2026 Giorgio Gonnella */

/*
   HyperLogLog: estimation of the number of distinct elements of a
   multiset, using a small fixed amount of memory: the estimator and its
   small-range correction (linear counting, for estimates up to 2.5 times
   the number of registers) of Flajolet et al. 2007; 64-bit hash values,
   as in Heule et al. 2013, so that no large-range correction is needed.

   The sketch is an array of 2^HLL_PRECISION one-byte registers (4 KiB);
   the relative standard error of the estimate is about
   1.04 / sqrt(2^HLL_PRECISION), i.e. 1.6%. Sketches of different parts of
   the data can be merged (e.g. computed by different threads).

   The elements are given as 64-bit hash values (see hll_hash()).
*/

#ifndef HYPERLOGLOG_H
#define HYPERLOGLOG_H

#include <stddef.h>
#include <stdint.h>

#define HLL_PRECISION 12
#define HLL_REGISTERS (1u << HLL_PRECISION)

typedef struct {
  uint8_t registers[HLL_REGISTERS];
} hll_t;

/* Initialize an empty sketch */
void hll_init(hll_t *self);

/* 64-bit hash of <len> bytes of <data> (suitable for hll_add()) */
uint64_t hll_hash(const void *data, size_t len);

/* Add an element, given its hash value */
void hll_add(hll_t *self, uint64_t hash);

/* Merge the sketch <other> into <self> */
void hll_merge(hll_t *self, const hll_t *other);

/* Estimated number of distinct elements */
double hll_estimate(const hll_t *self);

/* Unit tests */
void hll_test(void);

#endif
