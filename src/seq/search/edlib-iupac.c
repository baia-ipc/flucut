/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <pthread.h>
#include "edlib.h"
#include "edlib-iupac.h"
#include "iupac.h"

static EdlibEqualityPair pairs[64];
static pthread_once_t pairs_once = PTHREAD_ONCE_INIT;
static int n_pairs;

static void init_pairs(void) {
  assert(n_iupac_equalities <= sizeof(pairs) / sizeof(pairs[0]));
  for (size_t i = 0; i < n_iupac_equalities; i++) {
    pairs[i].first = iupac_equalities[i].first;
    pairs[i].second = iupac_equalities[i].second;
  }
  n_pairs = (int)n_iupac_equalities;
}

const EdlibEqualityPair* edlib_iupac_equalities(int *count) {
  pthread_once(&pairs_once, init_pairs);
  if (count != NULL)
    *count = n_pairs;
  return pairs;
}
