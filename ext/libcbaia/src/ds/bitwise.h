/* (c) 2017-2026 Giorgio Gonnella */

#ifndef BITWISE_H
#define BITWISE_H

#include <stdint.h>

/* Generic single bit setter/getter; the bit number can be up to the
   width of the type of E minus 1 */
#define GET_BIT(E, BITNUM) (((E) >> (BITNUM)) & 1)
#define SET_BIT(E, BITNUM) ((E) |= ((__typeof__(E))1 << (BITNUM)))
#define UNSET_BIT(E, BITNUM) ((E) &= ~((__typeof__(E))1 << (BITNUM)))

/* Mixes all bits of x, so that each bit of the result depends on all bits
   of x (finalizer of splitmix64, https://xoshiro.di.unimi.it/splitmix64.c);
   e.g. to hash an integer key, or to improve a weak hash value */
static inline uint64_t bitwise_mix64(uint64_t x) {
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

#endif
