/* (c) 2022-2026 Giorgio Gonnella */

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "alloc.h"
#include "bitwise.h"
#include "bloom-filter.h"
#include "test.h"

struct bloom_filter_t {
  unsigned char *bits;
  size_t mask;   /* number of bits - 1 */
};

bloom_filter_t* bloom_filter_new(size_t size_exp) {
  assert(size_exp + 3 < sizeof(size_t) * 8);
  bloom_filter_t *self;
  XMALLOC(self, sizeof(bloom_filter_t));
  self->mask = ((size_t)1 << (size_exp + 3)) - 1;
  XCALLOC(self->bits, (size_t)1 << size_exp, 1);
  return self;
}

void bloom_filter_delete(bloom_filter_t *self) {
  if (self == NULL)
    return;
  free(self->bits);
  free(self);
}

/* the second hash value: the halves of the first one swapped */
static inline uint64_t swap_halves(uint64_t hash) {
  return hash >> 32 | hash << 32;
}

static inline void set_bit(bloom_filter_t *self, uint64_t hash) {
  size_t addr = (size_t)hash & self->mask;
  self->bits[addr >> 3] |= (unsigned char)(1U << (addr & 7));
}

static inline bool get_bit(const bloom_filter_t *self, uint64_t hash) {
  size_t addr = (size_t)hash & self->mask;
  return (self->bits[addr >> 3] >> (addr & 7)) & 1;
}

void bloom_filter_insert(bloom_filter_t *self, size_t key) {
  assert(self != NULL);
  uint64_t hash = bitwise_mix64(key);
  set_bit(self, hash);
  set_bit(self, swap_halves(hash));
}

bool bloom_filter_query(const bloom_filter_t *self, size_t key) {
  assert(self != NULL);
  uint64_t hash = bitwise_mix64(key);
  return get_bit(self, hash) && get_bit(self, swap_halves(hash));
}

void bloom_filter_test(void) {
  /* the hash of the keys: the first output of splitmix64 with seed 0 */
  test_assert(bitwise_mix64(0x9e3779b97f4a7c15ULL) == 0xe220a8397b1dcdafULL);
  bloom_filter_t *bf = bloom_filter_new(16);
  for (size_t key = 0; key < 1000; key++)
    bloom_filter_insert(bf, key * 7);
  for (size_t key = 0; key < 1000; key++)
    test_assert(bloom_filter_query(bf, key * 7));
  /* 2^19 bits, 2000 bits set: few false positives expected */
  size_t false_positives = 0;
  for (size_t key = 0; key < 10000; key++)
    if (bloom_filter_query(bf, 1000000 + key))
      false_positives++;
  test_assert(false_positives < 100);
  bloom_filter_delete(bf);
}
