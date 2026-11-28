/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <pthread.h>
#include <stdbool.h>
#include <ctype.h>
#include "iupac.h"

/*
   Character pairs that match under IUPAC ambiguity rules: each IUPAC
   code with each of the bases it stands for, in both directions (a
   degenerate base can be in the primer, e.g. R, or in the read, e.g. N,
   which is also the padding character). The constant table is shared by
   all threads and by the pattern scanner.
*/
const iupac_pair_t iupac_equalities[] = {
  {'R', 'A'}, {'R', 'G'},
  {'Y', 'C'}, {'Y', 'T'},
  {'S', 'C'}, {'S', 'G'},
  {'W', 'A'}, {'W', 'T'},
  {'K', 'G'}, {'K', 'T'},
  {'M', 'A'}, {'M', 'C'},
  {'B', 'C'}, {'B', 'G'}, {'B', 'T'},
  {'D', 'A'}, {'D', 'G'}, {'D', 'T'},
  {'H', 'A'}, {'H', 'C'}, {'H', 'T'},
  {'V', 'A'}, {'V', 'C'}, {'V', 'G'},
  {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'},
  /* N (also the padding character) with the codes of the primer: a
     degenerate base of a primer overhanging the read is not an error */
  {'N', 'R'}, {'N', 'Y'}, {'N', 'S'}, {'N', 'W'}, {'N', 'K'}, {'N', 'M'},
  {'N', 'B'}, {'N', 'D'}, {'N', 'H'}, {'N', 'V'}, {'N', 'N'},
  {'U', 'T'},
};
const size_t n_iupac_equalities =
  sizeof(iupac_equalities) / sizeof(iupac_equalities[0]);

static iupac_equal_chars_t equal_chars[256];
static pthread_once_t equal_chars_once = PTHREAD_ONCE_INIT;

static void init_equal_chars(void) {
  for (int c = 0; c < 256; c++) {
    equal_chars[c].chars[0] = (unsigned char)c;
    equal_chars[c].n = 1;
  }
  for (size_t i = 0; i < n_iupac_equalities; i++) {
    unsigned char pair[2] = {(unsigned char)iupac_equalities[i].first,
                             (unsigned char)iupac_equalities[i].second};
    for (int side = 0; side < 2; side++) {
      iupac_equal_chars_t *a = equal_chars + pair[side];
      unsigned char b = pair[1 - side];
      bool known = false;
      for (int j = 0; j < a->n; j++)
        known = known || a->chars[j] == b;
      if (!known) {
        assert(a->n < IUPAC_MAX_EQUAL);
        a->chars[a->n++] = b;
      }
    }
  }
}

const iupac_equal_chars_t* iupac_equal_chars(void) {
  pthread_once(&equal_chars_once, init_equal_chars);
  return equal_chars;
}

static const struct { char code; const char *bases; } iupac_bases_table[] = {
  {'A', "A"}, {'C', "C"}, {'G', "G"}, {'T', "T"}, {'U', "T"},
  {'R', "AG"}, {'Y', "CT"}, {'S', "CG"}, {'W', "AT"}, {'K', "GT"},
  {'M', "AC"}, {'B', "CGT"}, {'D', "AGT"}, {'H', "ACT"},
  {'V', "ACG"}, {'N', "ACGT"},
};

const char* iupac_bases(char code) {
  code = (char)toupper((unsigned char)code);
  for (size_t i = 0; i < sizeof(iupac_bases_table) /
                       sizeof(iupac_bases_table[0]); i++)
    if (iupac_bases_table[i].code == code)
      return iupac_bases_table[i].bases;
  return NULL;
}

bool iupac_is_seq(const char *seq, size_t len) {
  assert(seq != NULL);
  if (len == 0)
    return false;
  for (size_t i = 0; i < len; i++)
    if (iupac_bases(seq[i]) == NULL)
      return false;
  return true;
}
