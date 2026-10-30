/* (c) 2026 Giorgio Gonnella */

/*
   Index of the k-mers of a set of sequences (e.g. reference sequences),
   for finding where a query sequence may align: the k-mers of the query
   which occur in the indexed sequences give candidate diagonals
   (position in the sequence minus position in the query), around which
   the query is then aligned by the caller (e.g. with edlib), instead of
   being aligned to the whole of every sequence.

   Implementation: the k-mers (of the bases A, C, G, T; k-mers containing
   another character are not indexed) are encoded with 2 bits per base;
   a table of 4^k + 1 offsets (4^k * 4 bytes: 4 MB for k = 10) points
   into the arrays of the occurrences (sequence and position), which are
   sorted by sequence and position. The index is built once and only read
   afterwards: it can be shared by several threads.

   An exact k-mer is needed for a candidate: a query with errors in every
   k-mer gives none.
*/

#ifndef KMER_INDEX_H
#define KMER_INDEX_H

#include <stddef.h>
#include <stdint.h>

/* the largest k (the table of the offsets has 4^k + 1 entries) */
#define KMER_INDEX_MAX_K 12

/* a candidate: a diagonal of an indexed sequence */
typedef struct {
  uint32_t seq;        /* index of the sequence (order of the input) */
  long diagonal;       /* position in the sequence of the first base of
                          the query, if the k-mer is part of an
                          ungapped alignment (can be negative) */
} kmer_index_hit_t;

typedef struct kmer_index_t kmer_index_t;

/*
   Index the k-mers of the <n> sequences <seqs> (lengths <lens>; upper or
   lower case); <k>: 1 to KMER_INDEX_MAX_K. The sequences are not needed
   after the construction. The total length must be less than 2^32.
*/
kmer_index_t* kmer_index_new(const char *const *seqs, const size_t *lens,
    size_t n, size_t k);

/* The k of the index */
size_t kmer_index_k(const kmer_index_t *self);

/* The number of occurrences of the k-mer <kmer> (k characters); 0 if it
   contains a character other than A, C, G, T */
size_t kmer_index_count(const kmer_index_t *self, const char *kmer);

/*
   The candidate diagonals of <query> (length <len>): for each k-mer of
   the query occurring at most <max_occ> times in the indexed sequences
   (0: no limit; a limit skips the k-mers of low complexity), the
   diagonals of its occurrences. They are stored in <*hits> (a buffer of
   <*cap> bytes, NULL and 0 at the first call, grown if needed with the
   XCHECKSPACE macros: it can be reused between calls and is freed by the
   caller), without duplicates, sorted by sequence and diagonal; returns
   their number.
*/
size_t kmer_index_diagonals(const kmer_index_t *self, const char *query,
    size_t len, size_t max_occ, kmer_index_hit_t **hits, size_t *cap);

/* Delete (NULL is accepted) */
void kmer_index_delete(kmer_index_t *self);

/* Unit tests */
void kmer_index_test(void);

#endif
