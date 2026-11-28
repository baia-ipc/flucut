/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "khash.h"
#include "alloc.h"
#include "file.h"
#include "cell-caller.h"
#include "strutils.h"
#include "test.h"

KHASH_MAP_INIT_INT64(counts, uint64_t)

struct cell_caller_t {
  khash_t(counts) *counts;
};

typedef struct {
  uint64_t barcode, count;
} barcode_count_t;

cell_caller_t* cell_caller_new(void) {
  cell_caller_t *self;
  XMALLOC(self, sizeof(cell_caller_t));
  self->counts = kh_init(counts);
  if (self->counts == NULL)
    fatal_printf("Failed allocating the counts of the barcodes\n");
  return self;
}

void cell_caller_delete(cell_caller_t *self) {
  if (self == NULL)
    return;
  kh_destroy(counts, self->counts);
  free(self);
}

static void add_count(cell_caller_t *self, uint64_t barcode,
    uint64_t count) {
  int ret;
  khiter_t k = kh_put(counts, self->counts, barcode, &ret);
  if (ret == -1)
    fatal_printf("Failed allocating the counts of the barcodes\n");
  if (ret != 0)
    kh_value(self->counts, k) = 0;
  kh_value(self->counts, k) += count;
}

void cell_caller_add(cell_caller_t *self, uint64_t barcode) {
  assert(self != NULL);
  add_count(self, barcode, 1);
}

void cell_caller_merge(cell_caller_t *self, cell_caller_t *other) {
  assert(self != NULL);
  assert(other != NULL);
  for (khiter_t k = kh_begin(other->counts); k != kh_end(other->counts);
       k++)
    if (kh_exist(other->counts, k))
      add_count(self, kh_key(other->counts, k),
                kh_value(other->counts, k));
  kh_clear(counts, other->counts);
}

/* decreasing count, then increasing barcode (for a stable output) */
static int cmp_counts(const void *a, const void *b) {
  const barcode_count_t *x = a, *y = b;
  if (x->count != y->count)
    return x->count < y->count ? 1 : -1;
  return x->barcode < y->barcode ? -1 : x->barcode > y->barcode;
}

/* the barcodes and their counts, sorted by cmp_counts, in a new array
   <bc> (to be freed by the caller); returns their number */
static size_t sorted_barcode_counts(const cell_caller_t *self,
    barcode_count_t **bc) {
  size_t n = kh_size(self->counts), i = 0;
  XMALLOC(*bc, sizeof(barcode_count_t) * (n + 1));
  for (khiter_t k = kh_begin(self->counts); k != kh_end(self->counts); k++)
    if (kh_exist(self->counts, k)) {
      (*bc)[i].barcode = kh_key(self->counts, k);
      (*bc)[i++].count = kh_value(self->counts, k);
    }
  qsort(*bc, n, sizeof(barcode_count_t), cmp_counts);
  return n;
}

size_t cell_caller_sorted_counts(const cell_caller_t *self,
    uint64_t **counts) {
  assert(self != NULL);
  assert(counts != NULL);
  barcode_count_t *bc;
  size_t n = sorted_barcode_counts(self, &bc);
  XMALLOC(*counts, sizeof(uint64_t) * (n + 1));
  for (size_t i = 0; i < n; i++)
    (*counts)[i] = bc[i].count;
  free(bc);
  return n;
}

whitelist_t* cell_caller_run(const cell_caller_t *self,
    size_t expected_cells, size_t length, const whitelist_t *known,
    const char *path, Logger *logger) {
  assert(self != NULL);
  assert(expected_cells > 0);
  barcode_count_t *bc;
  size_t n = sorted_barcode_counts(self, &bc), i;
  /* 99th percentile of the counts of the expected cells */
  size_t top = expected_cells < n ? expected_cells : n;
  size_t p99 = top == 0 ? 0 : (size_t)((double)(top - 1) * 0.01 + 0.5);
  /* rounded, as the cutoff of Cell Ranger 2 */
  uint64_t threshold = n == 0 ? 1 : (bc[p99].count + 5) / 10;
  if (threshold < 1)
    threshold = 1;
  size_t n_cells = 0;
  while (n_cells < n && bc[n_cells].count >= threshold)
    n_cells++;
  /* the known cells not selected */
  uint64_t *known_codes = NULL;
  size_t n_known = known == NULL ? 0 : whitelist_codes(known, &known_codes);
  uint64_t *codes;
  XMALLOC(codes, sizeof(uint64_t) * (n_cells + n_known + 1));
  uint64_t in_cells = 0, total = 0;
  for (i = 0; i < n; i++) {
    total += bc[i].count;
    if (i < n_cells) {
      codes[i] = bc[i].barcode;
      in_cells += bc[i].count;
    }
  }
  whitelist_t *selected = whitelist_new_from_codes(codes, n_cells, length);
  size_t n_added = 0, added_with_reads = 0;
  uint64_t in_added = 0;
  for (size_t k = 0; k < n_known; k++) {
    if (whitelist_contains(selected, known_codes[k]))
      continue;
    codes[n_cells + n_added++] = known_codes[k];
    khiter_t it = kh_get(counts, self->counts, known_codes[k]);
    if (it != kh_end(self->counts)) {
      added_with_reads++;
      in_added += kh_value(self->counts, it);
    }
  }
  whitelist_delete(selected);
  whitelist_t *cells = whitelist_new_from_codes(codes, n_cells + n_added,
                                                length);
  lprintf(logger, LOG_INFO, "Cell barcodes: %zu barcodes found exactly "
          "(%llu reads), %zu cells (at least %llu reads; %.1f%% of the "
          "reads)\n", n, (unsigned long long)total, n_cells,
          (unsigned long long)threshold,
          total > 0 ? 100.0 * (double)in_cells / (double)total : 0.0);
  if (known != NULL)
    lprintf(logger, LOG_INFO, "Known cells: %zu, of which %zu added to the "
            "cells (%zu of them found exactly, %.1f%% of the reads)\n",
            n_known, n_added, added_with_reads,
            total > 0 ? 100.0 * (double)in_added / (double)total : 0.0);
  if (path != NULL) {
    FILE *f = xfopen(path, NULL, "w", "cell barcodes", logger);
    char seq[WHITELIST_MAX_BARCODE_LENGTH + 1];
    for (i = 0; i < n_cells + n_added; i++) {
      uint64_t count = 0;
      if (i < n_cells)
        count = bc[i].count;
      else {
        khiter_t it = kh_get(counts, self->counts, codes[i]);
        if (it != kh_end(self->counts))
          count = kh_value(self->counts, it);
      }
      whitelist_decode(cells, codes[i], seq);
      fprintf(f, "%s\t%llu\n", seq, (unsigned long long)count);
    }
    xfclose(&f, logger);
  }
  free(known_codes);
  free(codes);
  free(bc);
  return cells;
}

/* --- tests --- */

void cell_caller_test(const char *tmpdir) {
  cell_caller_t *cc = cell_caller_new(), *cc2 = cell_caller_new();
  /* 3 cells with 100, 80 and 20 reads, background barcodes with 1 read */
  for (uint64_t b = 1; b <= 3; b++)
    for (size_t k = 0; k < (b == 1 ? 100 : b == 2 ? 80 : 20); k++)
      cell_caller_add(b % 2 ? cc : cc2, b);
  for (uint64_t b = 10; b < 50; b++)
    cell_caller_add(cc2, b);
  cell_caller_merge(cc, cc2);
  char *path = strcat_alloc(tmpdir, "/cbaia-tests-cells.txt");
  whitelist_t *cells = cell_caller_run(cc, 3, 4, NULL, path, NULL);
  test_eq_size(whitelist_size(cells), 3);
  test_assert(whitelist_contains(cells, 3));
  test_assert(!whitelist_contains(cells, 10));
  whitelist_delete(cells);
  FILE *f = fopen(path, "r");
  test_assert(f != NULL);
  char line[64];
  test_assert(fgets(line, sizeof(line), f) != NULL);
  test_eq_str(line, "AAAC\t100\n");
  fclose(f);
  unlink(path);
  free(path);
  uint64_t *sorted;
  test_eq_size(cell_caller_sorted_counts(cc, &sorted), 43);
  test_eq_size((size_t)sorted[0], 100);
  test_eq_size((size_t)sorted[42], 1);
  free(sorted);
  /* expected 1 cell: threshold 10 */
  cells = cell_caller_run(cc, 1, 4, NULL, NULL, NULL);
  test_eq_size(whitelist_size(cells), 3);
  whitelist_delete(cells);
  /* known cells: one already selected, one with 1 read, one without
     reads; all are added */
  uint64_t known_codes[] = {1, 10, 60};
  whitelist_t *known = whitelist_new_from_codes(known_codes, 3, 4);
  cells = cell_caller_run(cc, 3, 4, known, NULL, NULL);
  test_eq_size(whitelist_size(cells), 5);
  test_assert(whitelist_contains(cells, 10));
  test_assert(whitelist_contains(cells, 60));
  whitelist_delete(cells);
  whitelist_delete(known);
  cell_caller_delete(cc);
  cell_caller_delete(cc2);
  /* the threshold is rounded: 1.4 -> 1 (the barcodes with 1 read are
     selected), 1.5 -> 2 (they are not) */
  for (uint64_t top = 14; top <= 15; top++) {
    cc = cell_caller_new();
    for (size_t k = 0; k < top; k++)
      cell_caller_add(cc, 1);
    for (uint64_t b = 10; b < 20; b++)
      cell_caller_add(cc, b);
    cells = cell_caller_run(cc, 1, 4, NULL, NULL, NULL);
    test_eq_size(whitelist_size(cells), top == 14 ? 11 : 1);
    whitelist_delete(cells);
    cell_caller_delete(cc);
  }
}
