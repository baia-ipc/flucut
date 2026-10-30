/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <math.h>
#include <string.h>
#include "bitwise.h"
#include "hyperloglog.h"
#include "test.h"

void hll_init(hll_t *self) {
  assert(self != NULL);
  memset(self->registers, 0, HLL_REGISTERS);
}

uint64_t hll_hash(const void *data, size_t len) {
  assert(data != NULL || len == 0);
  /* FNV-1a, then mixed (the bits of FNV-1a are not uniform enough) */
  const unsigned char *p = data;
  uint64_t h = 0xcbf29ce484222325ULL;
  for (size_t i = 0; i < len; i++) {
    h ^= p[i];
    h *= 0x100000001b3ULL;
  }
  return bitwise_mix64(h);
}

void hll_add(hll_t *self, uint64_t hash) {
  assert(self != NULL);
  /* the first bits select the register; the rank is the position of the
     first 1 bit of the others */
  uint32_t reg = (uint32_t)(hash >> (64 - HLL_PRECISION));
  uint64_t rest = (hash << HLL_PRECISION) | (1ULL << (HLL_PRECISION - 1));
  uint8_t rank = (uint8_t)(__builtin_clzll(rest) + 1);
  if (rank > self->registers[reg])
    self->registers[reg] = rank;
}

void hll_merge(hll_t *self, const hll_t *other) {
  assert(self != NULL);
  assert(other != NULL);
  for (size_t i = 0; i < HLL_REGISTERS; i++)
    if (other->registers[i] > self->registers[i])
      self->registers[i] = other->registers[i];
}

double hll_estimate(const hll_t *self) {
  assert(self != NULL);
  const double m = HLL_REGISTERS;
  double sum = 0;
  size_t zeros = 0;
  for (size_t i = 0; i < HLL_REGISTERS; i++) {
    sum += ldexp(1.0, -(int)self->registers[i]);
    zeros += self->registers[i] == 0;
  }
  double alpha = 0.7213 / (1.0 + 1.079 / m);
  double estimate = alpha * m * m / sum;
  /* small range: linear counting (Flajolet et al. 2007); no large range
     correction with 64-bit hashes (Heule et al. 2013) */
  if (estimate <= 2.5 * m && zeros > 0)
    estimate = m * log(m / (double)zeros);
  return estimate;
}

/* --- unit tests --- */

void hll_test(void) {
  hll_t a, b;
  hll_init(&a);
  test_assert(hll_estimate(&a) == 0);
  /* small and large cardinalities: within 5% (3 standard errors) */
  size_t sizes[] = {100, 1000, 100000, 1000000};
  for (size_t s = 0; s < sizeof(sizes) / sizeof(sizes[0]); s++) {
    hll_init(&a);
    for (uint64_t i = 0; i < sizes[s]; i++) {
      hll_add(&a, hll_hash(&i, sizeof(i)));
      /* duplicates do not change the estimate */
      hll_add(&a, hll_hash(&i, sizeof(i)));
    }
    double e = hll_estimate(&a);
    test_assert(fabs(e - (double)sizes[s]) < 0.05 * (double)sizes[s]);
  }
  /* merge: union of two overlapping sets of 50000 elements (75000) */
  hll_init(&a);
  hll_init(&b);
  for (uint64_t i = 0; i < 50000; i++) {
    hll_add(&a, hll_hash(&i, sizeof(i)));
    uint64_t j = i + 25000;
    hll_add(&b, hll_hash(&j, sizeof(j)));
  }
  hll_merge(&a, &b);
  test_assert(fabs(hll_estimate(&a) - 75000.0) < 0.05 * 75000.0);
  /* the hash of strings */
  test_assert(hll_hash("ACGT", 4) == hll_hash("ACGT", 4));
  test_assert(hll_hash("ACGT", 4) != hll_hash("ACGA", 4));
}
