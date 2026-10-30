/* (c) 2022-2026 Giorgio Gonnella */

/*
   Bloom filter for integer keys (two hash functions, derived from the
   splitmix64 hash of the key).
*/

#ifndef BLOOM_FILTER_H
#define BLOOM_FILTER_H

#include <stdbool.h>
#include <stddef.h>

typedef struct bloom_filter_t bloom_filter_t;

/* Create an empty bloom filter of 2^<size_exp> bytes */
bloom_filter_t *bloom_filter_new(size_t size_exp);

/* Insert a key */
void bloom_filter_insert(bloom_filter_t *self, size_t key);

/* false if the key was certainly not inserted, true if it may have been */
bool bloom_filter_query(const bloom_filter_t *self, size_t key);

/* Delete a bloom filter (NULL is accepted) */
void bloom_filter_delete(bloom_filter_t *self);

/* Unit tests */
void bloom_filter_test(void);

#endif
