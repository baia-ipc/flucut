/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "kmer-index.h"
#include "test.h"

struct kmer_index_t {
  size_t k;
  uint32_t mask;
  /* the occurrences of the k-mer with code c are the entries offsets[c]
     to offsets[c + 1] - 1 of seq_of and pos_of */
  uint32_t *offsets;
  uint32_t *seq_of, *pos_of;
};

/* the 2-bit code of a base (either case), -1 for another character */
static int base_code(char c) {
  switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return -1;
  }
}

/*
   The k-mers of a sequence, one at a time: kmer_scan_next() returns true
   and sets <code> and <pos> for the next k-mer without characters other
   than A, C, G, T, false at the end of the sequence
*/
typedef struct {
  const char *seq;
  size_t len, i, valid, k;
  uint32_t code, mask;
} kmer_scan_t;

static void kmer_scan_init(kmer_scan_t *s, const char *seq, size_t len,
    size_t k, uint32_t mask) {
  s->seq = seq;
  s->len = len;
  s->i = 0;
  s->valid = 0;
  s->k = k;
  s->code = 0;
  s->mask = mask;
}

static bool kmer_scan_next(kmer_scan_t *s, uint32_t *code, size_t *pos) {
  while (s->i < s->len) {
    int b = base_code(s->seq[s->i++]);
    if (b < 0) {
      s->valid = 0;
      continue;
    }
    s->code = ((s->code << 2) | (uint32_t)b) & s->mask;
    if (++s->valid >= s->k) {
      *code = s->code;
      *pos = s->i - s->k;
      return true;
    }
  }
  return false;
}

kmer_index_t* kmer_index_new(const char *const *seqs, const size_t *lens,
    size_t n, size_t k) {
  assert(seqs != NULL || n == 0);
  assert(k >= 1 && k <= KMER_INDEX_MAX_K);
  kmer_index_t *self;
  XMALLOC(self, sizeof(kmer_index_t));
  self->k = k;
  size_t n_codes = (size_t)1 << (2 * k);
  self->mask = (uint32_t)(n_codes - 1);
  XCALLOC(self->offsets, n_codes + 1, sizeof(uint32_t));
  /* count the occurrences, then turn the counts into offsets */
  size_t total = 0;
  kmer_scan_t s;
  uint32_t code;
  size_t pos;
  for (size_t i = 0; i < n; i++) {
    assert(lens[i] < UINT32_MAX);
    kmer_scan_init(&s, seqs[i], lens[i], k, self->mask);
    while (kmer_scan_next(&s, &code, &pos)) {
      self->offsets[code + 1]++;
      total++;
    }
  }
  assert(total < UINT32_MAX);
  for (size_t c = 0; c < n_codes; c++)
    self->offsets[c + 1] += self->offsets[c];
  XMALLOC(self->seq_of, sizeof(uint32_t) * (total > 0 ? total : 1));
  XMALLOC(self->pos_of, sizeof(uint32_t) * (total > 0 ? total : 1));
  /* fill, in the order of the sequences and positions */
  uint32_t *next;
  XMALLOC(next, sizeof(uint32_t) * n_codes);
  memcpy(next, self->offsets, sizeof(uint32_t) * n_codes);
  for (size_t i = 0; i < n; i++) {
    kmer_scan_init(&s, seqs[i], lens[i], k, self->mask);
    while (kmer_scan_next(&s, &code, &pos)) {
      uint32_t e = next[code]++;
      self->seq_of[e] = (uint32_t)i;
      self->pos_of[e] = (uint32_t)pos;
    }
  }
  free(next);
  return self;
}

size_t kmer_index_k(const kmer_index_t *self) {
  assert(self != NULL);
  return self->k;
}

size_t kmer_index_count(const kmer_index_t *self, const char *kmer) {
  assert(self != NULL);
  assert(kmer != NULL);
  kmer_scan_t s;
  uint32_t code;
  size_t pos;
  kmer_scan_init(&s, kmer, self->k, self->k, self->mask);
  if (!kmer_scan_next(&s, &code, &pos))
    return 0;
  return self->offsets[code + 1] - self->offsets[code];
}

/* the order of the hits: by sequence, then diagonal */
static int compare_hits(const void *a, const void *b) {
  const kmer_index_hit_t *x = a, *y = b;
  if (x->seq != y->seq)
    return x->seq < y->seq ? -1 : 1;
  return (x->diagonal > y->diagonal) - (x->diagonal < y->diagonal);
}

size_t kmer_index_diagonals(const kmer_index_t *self, const char *query,
    size_t len, size_t max_occ, kmer_index_hit_t **hits, size_t *cap) {
  assert(self != NULL);
  assert(query != NULL || len == 0);
  assert(hits != NULL);
  assert(cap != NULL);
  size_t n = 0;
  kmer_scan_t s;
  uint32_t code;
  size_t pos;
  kmer_scan_init(&s, query, len, self->k, self->mask);
  while (kmer_scan_next(&s, &code, &pos)) {
    uint32_t lo = self->offsets[code], hi = self->offsets[code + 1];
    if (max_occ > 0 && hi - lo > max_occ)
      continue;
    XCHECKSPACEINCPERC(*hits, cap, sizeof(kmer_index_hit_t) * (n + hi - lo),
                       50);
    for (uint32_t e = lo; e < hi; e++) {
      (*hits)[n].seq = self->seq_of[e];
      (*hits)[n].diagonal = (long)self->pos_of[e] - (long)pos;
      n++;
    }
  }
  if (n == 0)
    return 0;
  qsort(*hits, n, sizeof(kmer_index_hit_t), compare_hits);
  size_t u = 1;
  for (size_t i = 1; i < n; i++)
    if (compare_hits(*hits + i, *hits + u - 1) != 0)
      (*hits)[u++] = (*hits)[i];
  return u;
}

void kmer_index_delete(kmer_index_t *self) {
  if (self == NULL)
    return;
  free(self->offsets);
  free(self->seq_of);
  free(self->pos_of);
  free(self);
}

void kmer_index_test(void) {
  /* two sequences; the k-mer ACGT occurs in both, and twice in the
     first; the N interrupts the k-mers */
  const char *seqs[] = {"ACGTTTACGTGG", "ccACGTNACG"};
  size_t lens[] = {12, 10};
  kmer_index_t *ix = kmer_index_new(seqs, lens, 2, 4);
  test_eq_size(kmer_index_k(ix), 4);
  test_eq_size(kmer_index_count(ix, "ACGT"), 3);
  test_eq_size(kmer_index_count(ix, "acgt"), 3);
  test_eq_size(kmer_index_count(ix, "TTTA"), 1);
  test_eq_size(kmer_index_count(ix, "GTNA"), 0);
  test_eq_size(kmer_index_count(ix, "CCCC"), 0);
  kmer_index_hit_t *hits = NULL;
  size_t cap = 0;
  /* the query TTACGT: TTAC at 4 of the first sequence (diagonal 4 - 0),
     TACG at 5 (5 - 1 = 4), ACGT at 0, 6 (diagonals -2, 4) and at 2 of
     the second (diagonal 0); the duplicates are removed */
  size_t n = kmer_index_diagonals(ix, "TTACGT", 6, 0, &hits, &cap);
  test_eq_size(n, 3);
  test_eq_size(hits[0].seq, 0);
  test_eq_long(hits[0].diagonal, -2);
  test_eq_size(hits[1].seq, 0);
  test_eq_long(hits[1].diagonal, 4);
  test_eq_size(hits[2].seq, 1);
  test_eq_long(hits[2].diagonal, 0);
  /* with at most 2 occurrences, ACGT (3) is skipped */
  n = kmer_index_diagonals(ix, "TTACGT", 6, 2, &hits, &cap);
  test_eq_size(n, 1);
  test_eq_long(hits[0].diagonal, 4);
  /* no k-mer, or a query shorter than k */
  test_eq_size(kmer_index_diagonals(ix, "GGGGGG", 6, 0, &hits, &cap), 0);
  test_eq_size(kmer_index_diagonals(ix, "ACG", 3, 0, &hits, &cap), 0);
  free(hits);
  kmer_index_delete(ix);
  /* an empty set of sequences */
  ix = kmer_index_new(NULL, NULL, 0, 3);
  test_eq_size(kmer_index_count(ix, "ACG"), 0);
  kmer_index_delete(ix);
}
