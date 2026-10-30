/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "lengths.h"
#include "test.h"

struct lengths_t {
  size_t n_groups;
  /* the counts: LENGTHS_N_BINS for each group */
  size_t *counts;
};

lengths_t* lengths_new(size_t n_groups) {
  assert(n_groups > 0);
  lengths_t *self;
  XMALLOC(self, sizeof(lengths_t));
  self->n_groups = n_groups;
  XCALLOC(self->counts, n_groups * LENGTHS_N_BINS, sizeof(size_t));
  return self;
}

size_t lengths_n_groups(const lengths_t *self) {
  assert(self != NULL);
  return self->n_groups;
}

/* the bin of a length */
static size_t bin_of(size_t len) {
  size_t bin = len / LENGTHS_BIN;
  return bin < LENGTHS_N_BINS ? bin : LENGTHS_N_BINS - 1;
}

void lengths_add(lengths_t *self, size_t group, size_t len) {
  assert(self != NULL);
  assert(group < self->n_groups);
  self->counts[group * LENGTHS_N_BINS + bin_of(len)]++;
}

void lengths_add_to_bin(lengths_t *self, size_t group, size_t bin,
    size_t n) {
  assert(self != NULL);
  assert(group < self->n_groups);
  assert(bin < LENGTHS_N_BINS);
  self->counts[group * LENGTHS_N_BINS + bin] += n;
}

void lengths_add_all(lengths_t *self, const lengths_t *other) {
  assert(self != NULL);
  assert(other != NULL);
  assert(self->n_groups == other->n_groups);
  for (size_t i = 0; i < self->n_groups * LENGTHS_N_BINS; i++)
    self->counts[i] += other->counts[i];
}

void lengths_clear(lengths_t *self) {
  assert(self != NULL);
  memset(self->counts, 0,
         sizeof(size_t) * self->n_groups * LENGTHS_N_BINS);
}

size_t lengths_count(const lengths_t *self, size_t group, size_t bin) {
  assert(self != NULL);
  assert(group < self->n_groups);
  assert(bin < LENGTHS_N_BINS);
  return self->counts[group * LENGTHS_N_BINS + bin];
}

void lengths_summarize(const lengths_t *self, size_t group,
    size_t expected, lengths_summary_t *summary) {
  assert(self != NULL);
  assert(group < self->n_groups);
  assert(summary != NULL);
  const size_t *c = self->counts + group * LENGTHS_N_BINS;
  memset(summary, 0, sizeof(lengths_summary_t));
  for (size_t b = 0; b < LENGTHS_N_BINS; b++) {
    summary->total += c[b];
    if (c[b] > summary->main_count) {
      summary->main_count = c[b];
      summary->main_bin = b;
    }
  }
  double limit = LENGTHS_SHORT_FRACTION * (double)expected;
  for (size_t b = 0; b < LENGTHS_N_BINS; b++) {
    if ((double)((b + 1) * LENGTHS_BIN) > limit)
      break;
    summary->short_count += c[b];
    if (c[b] > summary->short_bin_count) {
      summary->short_bin_count = c[b];
      summary->short_bin = b;
    }
  }
}

size_t lengths_peaks(const lengths_t *self, size_t group, size_t *bins,
    size_t max) {
  assert(self != NULL);
  assert(group < self->n_groups);
  assert(bins != NULL || max == 0);
  const size_t *c = self->counts + group * LENGTHS_N_BINS;
  size_t total = 0;
  for (size_t b = 0; b < LENGTHS_N_BINS; b++)
    total += c[b];
  size_t n = 0;
  for (size_t b = 0; b < LENGTHS_N_BINS; b++) {
    if (c[b] == 0 || (double)c[b] < LENGTHS_PEAK_FRACTION * (double)total)
      continue;
    size_t lo = b > LENGTHS_PEAK_RADIUS ? b - LENGTHS_PEAK_RADIUS : 0;
    size_t hi = b + LENGTHS_PEAK_RADIUS < LENGTHS_N_BINS ?
                b + LENGTHS_PEAK_RADIUS : LENGTHS_N_BINS - 1;
    bool is_peak = true;
    for (size_t o = lo; o <= hi && is_peak; o++)
      // (equal neighbors: the first one is the peak)
      if (c[o] > c[b] || (c[o] == c[b] && o < b))
        is_peak = false;
    if (!is_peak)
      continue;
    // insertion by decreasing count, keeping at most max
    size_t i = n < max ? n++ : max;
    while (i > 0 && c[bins[i - 1]] < c[b]) {
      if (i < max)
        bins[i] = bins[i - 1];
      i--;
    }
    if (i < max)
      bins[i] = b;
  }
  return n;
}

void lengths_write(const lengths_t *self, size_t group, const char *prefix,
    FILE *file) {
  assert(self != NULL);
  assert(group < self->n_groups);
  const size_t *c = self->counts + group * LENGTHS_N_BINS;
  for (size_t b = 0; b < LENGTHS_N_BINS; b++)
    if (c[b] > 0)
      fprintf(file, "%s%zu\t%zu\n", prefix, b * LENGTHS_BIN, c[b]);
}

void lengths_delete(lengths_t *self) {
  if (self == NULL)
    return;
  free(self->counts);
  free(self);
}

void lengths_test(void) {
  lengths_t *l = lengths_new(2);
  test_eq_size(lengths_n_groups(l), 2);
  /* group 0: 50 amplicons of 1750 bases, 10 of 690, 3 of 1700 (not
     short: above 0.9 of the expected 1754) */
  for (int i = 0; i < 50; i++)
    lengths_add(l, 0, 1754);
  for (int i = 0; i < 10; i++)
    lengths_add(l, 0, 693);
  for (int i = 0; i < 3; i++)
    lengths_add(l, 0, 1700);
  /* a length beyond the last bin */
  lengths_add(l, 1, LENGTHS_MAX + 500);
  test_eq_size(lengths_count(l, 0, 175), 50);
  test_eq_size(lengths_count(l, 0, 69), 10);
  test_eq_size(lengths_count(l, 1, LENGTHS_N_BINS - 1), 1);
  lengths_summary_t s;
  lengths_summarize(l, 0, 1754, &s);
  test_eq_size(s.total, 63);
  test_eq_size(s.main_bin, 175);
  test_eq_size(s.main_count, 50);
  test_eq_size(s.short_count, 10);
  test_eq_size(s.short_bin, 69);
  test_eq_size(s.short_bin_count, 10);
  /* the peaks: 1750 and 690 (1700 has 3 of 63 amplicons, above 2%,
     but 1750 is larger within 5 bins) */
  size_t peaks[4];
  test_eq_size(lengths_peaks(l, 0, peaks, 4), 2);
  test_eq_size(peaks[0], 175);
  test_eq_size(peaks[1], 69);
  test_eq_size(lengths_peaks(l, 0, peaks, 1), 1);
  test_eq_size(peaks[0], 175);
  /* the short amplicons are counted also when they are the majority */
  for (int i = 0; i < 100; i++)
    lengths_add(l, 0, 693);
  lengths_summarize(l, 0, 1754, &s);
  test_eq_size(s.main_bin, 69);
  test_eq_size(s.short_count, 110);
  /* the sum of two histograms, and an empty group */
  lengths_t *m = lengths_new(2);
  lengths_add_all(m, l);
  lengths_add_all(m, l);
  test_eq_size(lengths_count(m, 0, 175), 100);
  test_eq_size(lengths_count(m, 0, 69), 220);
  lengths_clear(m);
  lengths_summarize(m, 0, 1754, &s);
  test_eq_size(s.total, 0);
  test_eq_size(s.short_count, 0);
  lengths_delete(m);
  lengths_delete(l);
}
