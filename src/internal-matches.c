/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "internal-matches.h"
#include "logger.h"
#include "pattern-scan.h"
#include "test.h"

/* the greatest number of matches considered in a read, before the
   choice of those which do not overlap */
#define MAX_CANDIDATES 128

struct internal_matches_t {
  const scheme_t *scheme;
  /* for each primer, of the left end and then of the right end: the
     prepared pattern and the errors allowed */
  pattern_scan_t **scans;
  long *max_errors;
  /* the types of influenza whose primers are searched (see
     internal_matches_restrict_virus) */
  virus_set_t types;
};

internal_matches_t* internal_matches_new(const scheme_t *scheme,
    const long *left_max_errors, const long *right_max_errors) {
  assert(scheme != NULL);
  internal_matches_t *self;
  XMALLOC(self, sizeof(internal_matches_t));
  size_t n = scheme->n_left + scheme->n_right;
  self->scheme = scheme;
  self->types = VIRUS_SET_BOTH;
  XMALLOC(self->scans, sizeof(pattern_scan_t*) * n);
  XMALLOC(self->max_errors, sizeof(long) * n);
  for (size_t i = 0; i < n; i++) {
    bool left = i < scheme->n_left;
    const primer_t *p = left ? scheme->left + i :
                        scheme->right + (i - scheme->n_left);
    self->scans[i] = pattern_scan_new(p->seq, p->seq_len);
    self->max_errors[i] = left ? left_max_errors[i] :
                          right_max_errors[i - scheme->n_left];
  }
  return self;
}

void internal_matches_restrict_virus(internal_matches_t *self,
    virus_set_t types) {
  assert(self != NULL);
  self->types = types;
}

void internal_matches_delete(internal_matches_t *self) {
  if (self == NULL)
    return;
  size_t n = self->scheme->n_left + self->scheme->n_right;
  for (size_t i = 0; i < n; i++)
    pattern_scan_delete(self->scans[i]);
  free(self->scans);
  free(self->max_errors);
  free(self);
}

/* the order in which the matches are kept (see internal-matches.h) */
static int compare_quality(const void *a, const void *b) {
  const amplicon_match_t *x = a, *y = b;
  if (x->n_errors != y->n_errors)
    return x->n_errors < y->n_errors ? -1 : 1;
  if (x->at_left_end != y->at_left_end)
    return x->at_left_end ? -1 : 1;
  long lx = x->end - x->start, ly = y->end - y->start;
  if (lx != ly)
    return lx > ly ? -1 : 1;
  if (x->variant != y->variant)
    return x->variant < y->variant ? -1 : 1;
  return x->start < y->start ? -1 : (x->start > y->start ? 1 : 0);
}

static int compare_position(const void *a, const void *b) {
  long d = ((const amplicon_match_t*)a)->start -
           ((const amplicon_match_t*)b)->start;
  return d < 0 ? -1 : (d > 0 ? 1 : 0);
}

size_t internal_matches_find(internal_matches_t *self, const char *seq,
    long from, long to, amplicon_match_t *matches) {
  assert(self != NULL);
  assert(matches != NULL);
  if (to <= from)
    return 0;
  const scheme_t *sc = self->scheme;
  amplicon_match_t candidates[MAX_CANDIDATES];
  size_t n_candidates = 0, n = 0;
  for (size_t i = 0; i < sc->n_left + sc->n_right; i++) {
    bool left = i < sc->n_left;
    const primer_t *primer = left ? sc->left + i :
                             sc->right + (i - sc->n_left);
    if (!virus_set_has(self->types, primer->virus))
      continue;
    pattern_match_t found[PATTERN_SCAN_MAX_MATCHES];
    size_t n_found = pattern_scan_matches(self->scans[i], seq + from,
        (size_t)(to - from), self->max_errors[i], found);
    for (size_t j = 0; j < n_found && n_candidates < MAX_CANDIDATES; j++) {
      amplicon_match_t *m = candidates + n_candidates++;
      m->at_left_end = i < sc->n_left;
      m->variant = m->at_left_end ? i : i - sc->n_left;
      m->start = found[j].start + from;
      m->end = found[j].end + from;
      m->n_errors = found[j].n_errors;
    }
  }
  qsort(candidates, n_candidates, sizeof(amplicon_match_t),
        compare_quality);
  for (size_t i = 0; i < n_candidates && n < AMPLICON_MAX_INTERNAL; i++) {
    const amplicon_match_t *c = candidates + i;
    bool overlaps = false;
    for (size_t j = 0; j < n && !overlaps; j++)
      overlaps = c->start <= matches[j].end && matches[j].start <= c->end;
    if (!overlaps)
      matches[n++] = *c;
  }
  qsort(matches, n, sizeof(amplicon_match_t), compare_position);
  return n;
}

/* --- unit tests --- */

void internal_matches_test(void) {
  Logger *logger = new_logger(stderr, LOG_ERROR);
  scheme_t *sc = scheme_new_builtin(SCHEME_DEFAULT, logger);
  long left_errors[8], right_errors[8];
  assert(sc->n_left <= 8 && sc->n_right <= 8);
  for (size_t i = 0; i < sc->n_left; i++)
    left_errors[i] = 1;
  for (size_t i = 0; i < sc->n_right; i++)
    right_errors[i] = 2;
  internal_matches_t *im = internal_matches_new(sc, left_errors,
                                                right_errors);
  amplicon_match_t m[AMPLICON_MAX_INTERNAL];
  const char *filler = "GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAA";
  char read[512];

  /* a match with an error before an exact one: both are reported, in
     the order of the read (the better one does not hide the other) */
  snprintf(read, sizeof(read), "%sAGCAATAGCAGG%sAGCAAAAGCAGG%s", filler,
           filler, filler);
  test_eq_size(internal_matches_find(im, read, 0,
                                     (long)strlen(read), m), 2);
  test_eq_long(m[0].start, 32);
  test_eq_long(m[0].end, 43);
  test_eq_long(m[0].n_errors, 1);
  test_assert(m[0].at_left_end);
  test_eq_long(m[1].start, 76);
  test_eq_long(m[1].n_errors, 0);

  /* only the bases asked for are searched */
  test_eq_size(internal_matches_find(im, read, 44,
                                     (long)strlen(read), m), 1);
  test_eq_long(m[0].start, 76);
  test_eq_size(internal_matches_find(im, read, 33, 86, m), 0);
  test_eq_size(internal_matches_find(im, read, 50, 50, m), 0);

  /* a primer of the right end, as it is found at the end of an
     amplicon */
  snprintf(read, sizeof(read), "%sCCTTGTTTCTACTGATCACGCGT%s", filler,
           filler);
  test_eq_size(internal_matches_find(im, read, 0,
                                     (long)strlen(read), m), 1);
  test_assert(!m[0].at_left_end);
  test_eq_long(m[0].start, 32);
  test_eq_long(m[0].end, 54);

  /* more matches than are kept */
  read[0] = '\0';
  for (int i = 0; i < AMPLICON_MAX_INTERNAL + 2; i++) {
    strcat(read, "AGCAAAAGCAGG");
    strcat(read, "GGTTCCAAGGTT");
  }
  test_eq_size(internal_matches_find(im, read, 0,
                                     (long)strlen(read), m),
               AMPLICON_MAX_INTERNAL);

  internal_matches_delete(im);
  scheme_delete(sc);
  delete_logger(logger);
}
