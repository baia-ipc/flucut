/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "read-agreement.h"
#include "test.h"

void read_agreement_init(read_agreement_t *self, size_t n_refs) {
  assert(self != NULL);
  self->n_refs = n_refs;
  XCALLOC(self->agreeing, n_refs > 0 ? n_refs : 1, sizeof(size_t));
  memset(self->reads, 0, sizeof(self->reads));
}

void read_agreement_add(read_agreement_t *self, const int *segment,
    const size_t *read_votes, size_t n) {
  assert(self != NULL);
  if (n > READ_AGREEMENT_MAX_PROBES)
    n = READ_AGREEMENT_MAX_PROBES;
  /* the probes of a read may lie in several segments (a chimera): each
     segment is considered by itself */
  for (int s = 1; s <= SEGMENTS_N; s++) {
    size_t k = 0, first = 0;
    bool agree = true;
    for (size_t j = 0; j < n; j++) {
      if (segment[read_votes[j]] != s)
        continue;
      if (k == 0)
        first = read_votes[j];
      else if (read_votes[j] != first)
        agree = false;
      k++;
    }
    if (k < 2)
      continue;
    self->reads[s][k]++;
    if (agree)
      self->agreeing[first]++;
  }
}

size_t read_agreement_reads(const read_agreement_t *self, int s) {
  size_t n = 0;
  for (size_t k = 2; k <= READ_AGREEMENT_MAX_PROBES; k++)
    n += self->reads[s][k];
  return n;
}

double read_agreement_expected(const read_agreement_t *self, int s,
    double share) {
  double expected = 0, power = share;
  for (size_t k = 2; k <= READ_AGREEMENT_MAX_PROBES; k++) {
    power *= share;
    expected += (double)self->reads[s][k] * power;
  }
  return expected;
}

/* the votes of the segment <s> */
static size_t segment_votes(const int *segment, const size_t *votes,
    size_t n, int s) {
  size_t total = 0;
  for (size_t i = 0; i < n; i++)
    if (segment[i] == s)
      total += votes[i];
  return total;
}

void read_agreement_log(const read_agreement_t *self, const fasta_t *fa,
    const int *segment, const size_t *votes, const char *file,
    Logger *logger) {
  assert(self != NULL);
  assert(fa->n_records == self->n_refs);
  for (int s = 1; s <= SEGMENTS_N; s++) {
    size_t n_refs = 0, agreeing = 0;
    size_t total = segment_votes(segment, votes, self->n_refs, s);
    size_t n_reads = read_agreement_reads(self, s);
    double expected = 0;
    for (size_t i = 0; i < self->n_refs; i++) {
      if (segment[i] != s)
        continue;
      n_refs++;
      agreeing += self->agreeing[i];
      expected += read_agreement_expected(self, s,
          total > 0 ? (double)votes[i] / (double)total : 0);
    }
    if (n_refs < 2 || n_reads == 0)
      continue;
    char name[16];
    segments_name(s, name, sizeof(name));
    /* the references on which reads agree, with the reads expected */
    char list[1024] = "";
    size_t used = 0;
    for (size_t i = 0; i < self->n_refs; i++) {
      if (segment[i] != s || self->agreeing[i] == 0)
        continue;
      int w = snprintf(list + used, sizeof(list) - used,
          "%s%s %zu (%.1f)", used > 0 ? ", " : "", fa->records[i].name,
          self->agreeing[i], read_agreement_expected(self, s,
              (double)votes[i] / (double)total));
      if (w < 0 || (size_t)w >= sizeof(list) - used)
        break;
      used += (size_t)w;
    }
    lprintf(logger, LOG_INFO, "probes within a read of %s for %s: %zu "
            "reads with several voting probes, all for one reference in "
            "%zu (%.1f%%; %.1f%% expected if the probes voted "
            "independently of their read); by reference, with the reads "
            "expected: %s\n", name, file, n_reads, agreeing,
            100.0 * (double)agreeing / (double)n_reads,
            100.0 * expected / (double)n_reads, list);
  }
}

void read_agreement_free(read_agreement_t *self) {
  if (self == NULL)
    return;
  free(self->agreeing);
  self->agreeing = NULL;
}

void read_agreement_test(void) {
  /* segment 4 has the references 0 and 1, segment 6 the reference 2 */
  const int segment[] = {4, 4, 6};
  read_agreement_t a;
  read_agreement_init(&a, 3);
  const size_t same[] = {0, 0, 0};
  const size_t split[] = {0, 1, 0, 1};
  const size_t other[] = {1, 1};
  const size_t single[] = {1};
  /* a chimera: two probes agree in segment 4, one lies in segment 6 */
  const size_t chimera[] = {1, 2, 1};
  read_agreement_add(&a, segment, same, 3);
  read_agreement_add(&a, segment, split, 4);
  read_agreement_add(&a, segment, other, 2);
  read_agreement_add(&a, segment, single, 1);
  read_agreement_add(&a, segment, chimera, 3);
  test_eq_size(read_agreement_reads(&a, 4), 4);
  test_eq_size(read_agreement_reads(&a, 6), 0);
  test_eq_size(a.agreeing[0], 1);
  test_eq_size(a.agreeing[1], 2);
  test_eq_size(a.agreeing[2], 0);
  test_eq_size(a.reads[4][2], 2);
  test_eq_size(a.reads[4][3], 1);
  test_eq_size(a.reads[4][4], 1);
  /* expected for a share of 0.5: 2 * 0.25 + 0.125 + 0.0625 */
  test_eq(read_agreement_expected(&a, 4, 0.5), 0.6875, double, "%g");
  test_eq(read_agreement_expected(&a, 4, 1.0), 4.0, double, "%g");
  test_eq(read_agreement_expected(&a, 4, 0.0), 0.0, double, "%g");
  read_agreement_free(&a);
}
