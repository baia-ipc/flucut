/* (c) 2026 Giorgio Gonnella */

#ifndef EDLIB_IUPAC_H
#define EDLIB_IUPAC_H

#include "edlib.h"

/* Edlib adapter for libcbaia's dependency-free IUPAC equivalence table. */
const EdlibEqualityPair* edlib_iupac_equalities(int *n_pairs);

#endif
