/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "logger.h"
#include "pattern-scan.h"
#include "primer-search.h"
#include "primers.h"
#include "terminus-matches.h"
#include "test.h"

/* the greatest number of matches of the termini considered in a read,
   before the choice of those which do not overlap */
#define MAX_CANDIDATES 128

/* a whole primer (with its tail), as it is found on the strand of the
   read */
typedef struct {
  char *seq;
  size_t len;
  long max_errors;
} oligo_t;

/* a terminus searched, as it is found on the strand of the read */
typedef struct {
  char *seq;
  size_t len;
  long max_errors;
  bool start;          /* the start of an amplicon (else the end of one) */
  pattern_scan_t *scan;
  /* the primers with a tail which have this terminus: it is their
     suffix at the start of an amplicon, their prefix at the end of one */
  oligo_t *primers;
  size_t n_primers;
  /* the types of influenza whose primers have this terminus */
  virus_set_t types;
} pattern_t;

struct terminus_matches_t {
  pattern_t *patterns;
  size_t n_patterns;
  size_t max_primers;  /* room in the primers of each pattern */
  /* the types of influenza whose termini are searched (see
     terminus_matches_restrict_virus) */
  virus_set_t types;
  primer_search_t *search;
};

/* the errors allowed in a sequence of length <len>: a number of bases,
   or a fraction of its length if smaller than 1 */
static long max_errors_of(double max_errors, size_t len) {
  if (max_errors >= 1.0)
    return (long)max_errors;
  long n = (long)(max_errors * (double)len);
  return n > 0 ? n : 0;
}

static char* copy_seq(const char *seq, size_t len) {
  char *copy;
  XMALLOC(copy, len + 1);
  memcpy(copy, seq, len);
  copy[len] = '\0';
  return copy;
}

/* the pattern <seq> (length <len>), added unless it is already there */
static pattern_t* add_pattern(terminus_matches_t *self, const char *seq,
    size_t len, bool start, double max_errors) {
  for (size_t i = 0; i < self->n_patterns; i++) {
    pattern_t *p = self->patterns + i;
    if (p->start == start && p->len == len && memcmp(p->seq, seq, len) == 0)
      return p;
  }
  pattern_t *p = self->patterns + self->n_patterns++;
  p->seq = copy_seq(seq, len);
  p->len = len;
  p->max_errors = max_errors_of(max_errors, len);
  if (p->max_errors > TERMINUS_MATCHES_MAX_ERRORS)
    p->max_errors = TERMINUS_MATCHES_MAX_ERRORS;
  if (p->max_errors >= (long)len)
    p->max_errors = (long)len - 1;
  p->start = start;
  p->scan = pattern_scan_new(seq, len);
  XMALLOC(p->primers, sizeof(oligo_t) * self->max_primers);
  p->n_primers = 0;
  p->types = 0;   /* (the primers added to it set theirs) */
  return p;
}

/* add the primer <seq> (length <len>) to those of the terminus <p>,
   unless it has no tail or is already there */
static void add_primer(pattern_t *p, const char *seq, size_t len,
    double max_errors) {
  if (len <= p->len)
    return;
  for (size_t i = 0; i < p->n_primers; i++)
    if (p->primers[i].len == len &&
        memcmp(p->primers[i].seq, seq, len) == 0)
      return;
  oligo_t *o = p->primers + p->n_primers++;
  o->seq = copy_seq(seq, len);
  o->len = len;
  o->max_errors = max_errors_of(max_errors, len);
}

/* add the terminus of the primer <p> as it is in the scheme (<start>:
   whether it is then the start of an amplicon) and its reverse
   complement, which is the opposite, each with the primer on its
   strand */
static void add_terminus(terminus_matches_t *self, const primer_t *p,
    bool start, double max_errors) {
  const char *terminus = p->seq + p->terminus_from;
  char *rc, *rc_terminus;
  XMALLOC(rc, p->seq_len + 1);
  XMALLOC(rc_terminus, p->terminus_len + 1);
  primers_reverse_complement(p->seq, p->seq_len, rc);
  primers_reverse_complement(terminus, p->terminus_len, rc_terminus);
  pattern_t *forward = add_pattern(self, terminus, p->terminus_len, start,
                                   max_errors);
  add_primer(forward, p->seq, p->seq_len, max_errors);
  forward->types |= VIRUS_SET_OF(p->virus);
  pattern_t *reverse = add_pattern(self, rc_terminus, p->terminus_len,
                                   !start, max_errors);
  add_primer(reverse, rc, p->seq_len, max_errors);
  reverse->types |= VIRUS_SET_OF(p->virus);
  free(rc);
  free(rc_terminus);
}

terminus_matches_t* terminus_matches_new(const scheme_t *scheme,
    double max_errors) {
  assert(scheme != NULL);
  terminus_matches_t *self;
  XMALLOC(self, sizeof(terminus_matches_t));
  self->max_primers = scheme->n_left + scheme->n_right;
  XMALLOC(self->patterns, sizeof(pattern_t) * 2 * self->max_primers);
  self->n_patterns = 0;
  self->types = VIRUS_SET_BOTH;
  /* (the primers of the right end are stored as they are found at the
     end of an amplicon, see scheme.h) */
  for (size_t i = 0; i < scheme->n_left; i++)
    add_terminus(self, scheme->left + i, true, max_errors);
  for (size_t i = 0; i < scheme->n_right; i++)
    add_terminus(self, scheme->right + i, false, max_errors);
  self->search = primer_search_new();
  return self;
}

void terminus_matches_restrict_virus(terminus_matches_t *self,
    virus_set_t types) {
  assert(self != NULL);
  self->types = types;
}

void terminus_matches_delete(terminus_matches_t *self) {
  if (self == NULL)
    return;
  for (size_t i = 0; i < self->n_patterns; i++) {
    pattern_t *p = self->patterns + i;
    for (size_t j = 0; j < p->n_primers; j++)
      free(p->primers[j].seq);
    free(p->primers);
    pattern_scan_delete(p->scan);
    free(p->seq);
  }
  free(self->patterns);
  primer_search_delete(self->search);
  free(self);
}

/* by errors, then by pattern, then by position */
static int compare_quality(const void *a, const void *b) {
  const terminus_match_t *x = a, *y = b;
  if (x->n_errors != y->n_errors)
    return x->n_errors < y->n_errors ? -1 : 1;
  if (x->pattern != y->pattern)
    return x->pattern < y->pattern ? -1 : 1;
  return x->from < y->from ? -1 : (x->from > y->from ? 1 : 0);
}

static int compare_position(const void *a, const void *b) {
  const terminus_match_t *x = a, *y = b;
  return x->from < y->from ? -1 : (x->from > y->from ? 1 : 0);
}

/* all the matches of the pattern <index> in <seq> (length <len>),
   appended to <matches> (<n> of them, room for MAX_CANDIDATES) */
static void pattern_matches(terminus_matches_t *self, size_t index,
    const char *seq, size_t len, terminus_match_t *matches, size_t *n) {
  const pattern_t *p = self->patterns + index;
  pattern_match_t found[PATTERN_SCAN_MAX_MATCHES];
  size_t n_found = pattern_scan_matches(p->scan, seq, len, p->max_errors,
                                        found);
  for (size_t i = 0; i < n_found && *n < MAX_CANDIDATES; i++) {
    terminus_match_t *m = matches + (*n)++;
    m->from = found[i].start;
    m->to = found[i].end;
    m->n_errors = found[i].n_errors;
    m->pattern = index;
    m->start = p->start;
  }
}

size_t terminus_matches_find(terminus_matches_t *self, const char *seq,
    size_t len, terminus_match_t *matches) {
  assert(self != NULL);
  assert(matches != NULL);
  terminus_match_t candidates[MAX_CANDIDATES];
  size_t n_candidates = 0, n = 0;
  for (size_t i = 0; i < self->n_patterns; i++)
    if ((self->patterns[i].types & self->types) != 0)
      pattern_matches(self, i, seq, len, candidates, &n_candidates);
  qsort(candidates, n_candidates, sizeof(terminus_match_t), compare_quality);
  for (size_t i = 0; i < n_candidates && n < TERMINUS_MATCHES_MAX; i++) {
    const terminus_match_t *c = candidates + i;
    bool overlaps = false;
    for (size_t j = 0; j < n && !overlaps; j++)
      overlaps = c->from <= matches[j].to && matches[j].from <= c->to;
    if (!overlaps)
      matches[n++] = *c;
  }
  qsort(matches, n, sizeof(terminus_match_t), compare_position);
  return n;
}

long terminus_matches_primer_bound(terminus_matches_t *self,
    const char *seq, size_t len, const terminus_match_t *m) {
  const pattern_t *p = self->patterns + m->pattern;
  long bound = m->start ? m->from : m->to;
  for (size_t i = 0; i < p->n_primers; i++) {
    const oligo_t *o = p->primers + i;
    long reach = (long)(o->len - p->len) + o->max_errors;
    long lo = m->start ? m->from - reach : m->from;
    long hi = m->start ? m->to + 1 : m->to + 1 + reach;
    if (lo < 0)
      lo = 0;
    if (hi > (long)len)
      hi = (long)len;
    primer_hit_t hit;
    if (!primer_search_run_seq(self->search, o->seq, o->len, seq + lo,
                               (size_t)(hi - lo), o->max_errors, false,
                               false, &hit))
      continue;
    if (m->start && hit.start + lo < bound)
      bound = hit.start + lo;
    if (!m->start && hit.end + lo > bound)
      bound = hit.end + lo;
  }
  return bound;
}


/* --- unit tests --- */

void terminus_matches_test(void) {
  Logger *logger = new_logger(stderr, LOG_ERROR);
  scheme_t *sc = scheme_new_builtin(SCHEME_DEFAULT, logger);
  /* (30%: 3 errors in a terminus, limited to 2) */
  terminus_matches_t *tm = terminus_matches_new(sc, 0.3);
  terminus_match_t m[TERMINUS_MATCHES_MAX];
  const char *filler = "GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAA";
  char read[512];

  /* the four forms, in the order of the read: the start of an amplicon
     (Uni12), the end of one (the reverse complement of Uni13), then the
     same on the other strand */
  snprintf(read, sizeof(read), "%s%s%s%s%s%s%s%s", "AGCGAAAGCAGG", filler,
           "CCTTGTTTCTACT", filler, "AGTAGAAACAAGG", filler,
           "CCTGCTTTTGCT", filler);
  test_eq_size(terminus_matches_find(tm, read, strlen(read), m), 4);
  test_eq_long(m[0].from, 0);
  test_eq_long(m[0].to, 11);
  test_assert(m[0].start);
  test_eq_long(m[1].from, 44);
  test_eq_long(m[1].to, 56);
  test_assert(!m[1].start);
  test_eq_long(m[2].from, 89);
  test_assert(m[2].start);
  test_eq_long(m[3].from, 134);
  test_assert(!m[3].start);
  test_eq_long(m[3].n_errors, 0);

  /* two errors are found, three are not, whatever is allowed in a
     primer; a better match after a worse one does not hide it */
  snprintf(read, sizeof(read), "%s%s%s%s%s%s", filler, "AGCTAATGCAGG",
           filler, "AGCTATTGCAGG", filler, "AGCAAAAGCAGG");
  test_eq_size(terminus_matches_find(tm, read, strlen(read), m), 2);
  test_eq_long(m[0].from, 32);
  test_eq_long(m[0].n_errors, 2);
  test_eq_long(m[1].from, 120);
  test_eq_long(m[1].n_errors, 0);

  /* the tail of a primer is on the outer side of its terminus: before
     the start of an amplicon, after the end of one; without the tail in
     the read, or for a primer without one, the bound of the match */
  snprintf(read, sizeof(read), "%s%s%s%s%s", filler,
           "ACGCGTGATCAGTAGAAACAAGG", filler, "CCTTGTTTCTACTGATCACGCGT",
           filler);
  test_eq_size(terminus_matches_find(tm, read, strlen(read), m), 2);
  test_eq_long(m[0].from, 42);
  test_eq_long(terminus_matches_primer_bound(tm, read, strlen(read), m),
               32);
  test_eq_long(m[1].to, 99);
  test_eq_long(terminus_matches_primer_bound(tm, read, strlen(read),
                                             m + 1), 109);
  snprintf(read, sizeof(read), "%s%s%s%s", filler, "AGTAGAAACAAGG", filler,
           "AGCAAAAGCAGG");
  test_eq_size(terminus_matches_find(tm, read, strlen(read), m), 2);
  test_eq_long(terminus_matches_primer_bound(tm, read, strlen(read), m),
               32);
  test_eq_long(terminus_matches_primer_bound(tm, read, strlen(read),
                                             m + 1), 77);

  test_eq_size(terminus_matches_find(tm, "", 0, m), 0);

  /* Restricted to one type, the termini of the other one are not
     searched. The search is exact here (0 errors): the conserved
     termini of the two types are so alike that with the 2 errors
     allowed above the terminus of influenza A also matches the one of
     influenza B (AGCAAAAGCAGG against AGCAGAAGCAGAGC: two
     differences), which is why the primers of a type are searched
     alone once the references say which type a sample is. */
  terminus_matches_t *exact = terminus_matches_new(sc, 0.0);
  snprintf(read, sizeof(read), "%s%s%s%s%s", filler, "AGCAAAAGCAGG",
           filler, "AGCAGAAGCAGAGC", filler);
  test_eq_size(terminus_matches_find(exact, read, strlen(read), m), 2);
  terminus_matches_restrict_virus(exact, VIRUS_SET_A);
  test_eq_size(terminus_matches_find(exact, read, strlen(read), m), 1);
  test_eq_long(m[0].from, 32);
  terminus_matches_restrict_virus(exact, VIRUS_SET_B);
  test_eq_size(terminus_matches_find(exact, read, strlen(read), m), 1);
  test_eq_long(m[0].from, 76);
  terminus_matches_restrict_virus(exact, VIRUS_SET_BOTH);
  test_eq_size(terminus_matches_find(exact, read, strlen(read), m), 2);
  terminus_matches_delete(exact);
  terminus_matches_delete(tm);
  scheme_delete(sc);
  delete_logger(logger);
}
