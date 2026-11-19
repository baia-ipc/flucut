/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "chimera.h"
#include "fasta.h"
#include "primer-search.h"
#include "primers.h"
#include "test.h"

/* the greatest number of matches of the termini considered in a read */
#define MAX_MATCHES 16

/* a terminus searched, as it is found on the strand of the read */
typedef struct {
  char *seq;
  size_t len;
  long max_errors;
  bool start;          /* the start of an amplicon (else the end of one) */
} pattern_t;

/* a match of a terminus in a read (0-based, inclusive) */
typedef struct {
  long from, to;
  bool start;
} match_t;

struct chimera_finder_t {
  pattern_t *patterns;
  size_t n_patterns;
  long window;
  primer_search_t *search;
  reference_set_t *refs;
};

/* the errors allowed in a terminus (see chimera_finder_new()) */
static long max_errors_of(double max_errors, size_t len) {
  if (max_errors >= 1.0)
    return (long)max_errors;
  long n = (long)(max_errors * (double)len);
  return n > 0 ? n : 0;
}

/* add the pattern <seq> (length <len>), unless it is already there */
static void add_pattern(chimera_finder_t *self, const char *seq,
    size_t len, bool start, double max_errors) {
  for (size_t i = 0; i < self->n_patterns; i++) {
    const pattern_t *p = self->patterns + i;
    if (p->start == start && p->len == len && memcmp(p->seq, seq, len) == 0)
      return;
  }
  pattern_t *p = self->patterns + self->n_patterns++;
  XMALLOC(p->seq, len + 1);
  memcpy(p->seq, seq, len);
  p->seq[len] = '\0';
  p->len = len;
  p->max_errors = max_errors_of(max_errors, len);
  p->start = start;
}

/* add the terminus of the primer <p> as it is in the scheme (<start>:
   whether it is then the start of an amplicon) and its reverse
   complement, which is the opposite */
static void add_terminus(chimera_finder_t *self, const primer_t *p,
    bool start, double max_errors) {
  const char *terminus = p->seq + p->terminus_from;
  char *rc;
  XMALLOC(rc, p->terminus_len + 1);
  primers_reverse_complement(terminus, p->terminus_len, rc);
  add_pattern(self, terminus, p->terminus_len, start, max_errors);
  add_pattern(self, rc, p->terminus_len, !start, max_errors);
  free(rc);
}

chimera_finder_t* chimera_finder_new(const scheme_t *scheme,
    double max_errors, long window) {
  assert(scheme != NULL);
  assert(window >= 0);
  chimera_finder_t *self;
  XMALLOC(self, sizeof(chimera_finder_t));
  XMALLOC(self->patterns,
          sizeof(pattern_t) * 2 * (scheme->n_left + scheme->n_right));
  self->n_patterns = 0;
  /* (the primers of the right end are stored as they are found at the
     end of an amplicon, see scheme.h) */
  for (size_t i = 0; i < scheme->n_left; i++)
    add_terminus(self, scheme->left + i, true, max_errors);
  for (size_t i = 0; i < scheme->n_right; i++)
    add_terminus(self, scheme->right + i, false, max_errors);
  self->window = window;
  self->search = primer_search_new();
  self->refs = NULL;
  return self;
}

void chimera_finder_use_references(chimera_finder_t *self,
    reference_set_t *refs) {
  assert(self != NULL);
  assert(refs != NULL);
  self->refs = refs;
}

void chimera_finder_delete(chimera_finder_t *self) {
  if (self == NULL)
    return;
  for (size_t i = 0; i < self->n_patterns; i++)
    free(self->patterns[i].seq);
  free(self->patterns);
  primer_search_delete(self->search);
  free(self);
}

/* the best match of a terminus in <seq> from <lo> to <hi> - 1 (fewest
   errors, then the first pattern), into <m>; false if there is none */
static bool best_match(chimera_finder_t *self, const char *seq, long lo,
    long hi, match_t *m) {
  bool found = false;
  long best_errors = 0;
  for (size_t i = 0; i < self->n_patterns; i++) {
    const pattern_t *p = self->patterns + i;
    primer_hit_t hit;
    if (!primer_search_run_seq(self->search, p->seq, p->len, seq + lo,
                               (size_t)(hi - lo), p->max_errors, false,
                               false, &hit))
      continue;
    if (found && hit.n_errors >= best_errors)
      continue;
    found = true;
    best_errors = hit.n_errors;
    m->from = hit.start + lo;
    m->to = hit.end + lo;
    m->start = p->start;
  }
  return found;
}

static int compare_matches(const void *a, const void *b) {
  long d = ((const match_t*)a)->from - ((const match_t*)b)->from;
  return d < 0 ? -1 : (d > 0 ? 1 : 0);
}

/* the matches of the termini in <seq> (length <len>), which do not
   overlap, sorted by position, into <matches> (room for MAX_MATCHES);
   their number is returned */
static size_t find_matches(chimera_finder_t *self, const char *seq,
    size_t len, match_t *matches) {
  /* the intervals still to search: the best match of an interval is
     taken, then both sides of it are searched */
  long stack[2 * MAX_MATCHES + 2][2];
  size_t n_stack = 0, n = 0;
  stack[n_stack][0] = 0;
  stack[n_stack++][1] = (long)len;
  while (n_stack > 0 && n < MAX_MATCHES) {
    long lo = stack[--n_stack][0], hi = stack[n_stack][1];
    match_t m;
    if (hi <= lo || !best_match(self, seq, lo, hi, &m))
      continue;
    matches[n++] = m;
    stack[n_stack][0] = lo;
    stack[n_stack++][1] = m.from;
    stack[n_stack][0] = m.to + 1;
    stack[n_stack++][1] = hi;
  }
  qsort(matches, n, sizeof(match_t), compare_matches);
  return n;
}

/* locate in the references the bases of <seq> (length <len>) from
   <from>, at most a probe of them, limited to the read */
static void locate(chimera_finder_t *self, const char *seq, size_t len,
    long from, long to, reference_hit_t *hit) {
  hit->mapped = false;
  if (from < 0)
    from = 0;
  if (to > (long)len)
    to = (long)len;
  if (to - from >= REFERENCE_MIN_PROBE_LEN)
    reference_set_locate(self->refs, seq + from, (size_t)(to - from), hit);
}

/*
   Whether the bases of <seq> from <a_from> to <a_to> - 1 and those from
   <b_from> to <b_to> - 1 (after the first ones) are located on the same
   reference and strand, at the distance they have in the read, give or
   take <tolerance>; <b> is the location of the second ones.
*/
static bool is_continuous(chimera_finder_t *self, const char *seq,
    size_t len, long a_from, long a_to, long b_from, long b_to,
    long tolerance) {
  reference_hit_t a, b;
  locate(self, seq, len, a_from, a_to, &a);
  locate(self, seq, len, b_from, b_to, &b);
  if (!a.mapped || !b.mapped || a.ref_index != b.ref_index ||
      a.reverse != b.reverse)
    return false;
  long in_read = b_from - (a_from < 0 ? 0 : a_from);
  long in_ref = b.reverse ? a.ref_start - b.ref_start :
                b.ref_start - a.ref_start;
  long diff = in_ref - in_read;
  return diff <= tolerance && diff >= -tolerance;
}

/*
   Whether the read stops at the match <m> (see chimera.h): the bases on
   the amplicon side of the match are located in the references, and
   those on the other side are not found next to them. The other side is
   given CHIMERA_OUTER_PROBES chances, the probe next to the match and
   those after it, as the segment can differ there from the references
   (e.g. at the cleavage site of the HA, just before its Uni12-like
   site): a read which runs through a site of the segment must not be
   cut there.
*/
static bool is_discontinuous(chimera_finder_t *self, const char *seq,
    size_t len, const match_t *m) {
  long probe = (long)reference_set_probe_len(self->refs);
  long tolerance = CHIMERA_CONTINUITY_TOLERANCE;
  long before = m->from, after = m->to + 1;
  reference_hit_t inner;
  if (m->start)
    locate(self, seq, len, after, after + probe, &inner);
  else
    locate(self, seq, len, before - probe, before, &inner);
  if (!inner.mapped)
    return false;
  for (long i = 0; i < CHIMERA_OUTER_PROBES; i++) {
    /* (the probe of the other side, <i> probes away from the match) */
    long a_from = before - probe, b_from = after;
    if (m->start)
      a_from -= i * probe;
    else
      b_from += i * probe;
    if (is_continuous(self, seq, len, a_from, a_from + probe, b_from,
                      b_from + probe, (i + 1) * tolerance))
      return false;
  }
  return true;
}

/* end the current part at <to> and begin the next one at <from>; false
   if there is no room for another part */
static bool cut(chimera_part_t *parts, size_t *n, long to, long from) {
  if (*n + 1 >= CHIMERA_MAX_PARTS)
    return false;
  chimera_part_t *last = parts + *n - 1;
  /* (a part without any base, e.g. before a match at the start of the
     read, is not kept) */
  if (to < last->from)
    last->from = from;
  else {
    last->to = to;
    parts[(*n)++].from = from;
  }
  return true;
}

size_t chimera_finder_run(chimera_finder_t *self, const char *seq,
    size_t len, chimera_part_t *parts) {
  assert(self != NULL);
  assert(seq != NULL || len == 0);
  assert(parts != NULL);
  if (len == 0)
    return 0;
  match_t matches[MAX_MATCHES];
  size_t n_matches = find_matches(self, seq, len, matches);
  size_t n = 1;
  parts[0].from = 0;
  for (size_t i = 0; i < n_matches; i++) {
    const match_t *m = matches + i;
    bool ok = true;
    if (!m->start && i + 1 < n_matches && matches[i + 1].start &&
        matches[i + 1].from - m->to - 1 <= CHIMERA_MAX_GAP) {
      /* a pair: the gap is dropped */
      ok = cut(parts, &n, m->to, matches[i + 1].from);
      i++;
    } else if (self->refs != NULL) {
      /* a single match, not the primer of an end of the read */
      bool inside = m->start ? m->from > self->window :
                    m->to < (long)len - 1 - self->window;
      if (inside && is_discontinuous(self, seq, len, m))
        ok = m->start ? cut(parts, &n, m->from - 1, m->from) :
             cut(parts, &n, m->to, m->to + 1);
    }
    if (!ok)
      break;
  }
  parts[n - 1].to = (long)len - 1;
  if (parts[n - 1].to < parts[n - 1].from)
    n--;
  return n;
}

/* --- unit tests --- */

/* a deterministic pseudo-random sequence of <len> bases into <seq> */
static void random_seq(char *seq, size_t len, unsigned long *state) {
  for (size_t i = 0; i < len; i++) {
    *state = *state * 6364136223846793005UL + 1442695040888963407UL;
    seq[i] = "ACGT"[(*state >> 33) & 3];
  }
  seq[len] = '\0';
}

/* split <read> and check the parts (<n> pairs of from and to) */
static void check_parts(chimera_finder_t *f, const char *read, size_t n,
    const long *expected) {
  chimera_part_t parts[CHIMERA_MAX_PARTS];
  test_eq_size(chimera_finder_run(f, read, strlen(read), parts), n);
  for (size_t i = 0; i < n; i++) {
    test_eq_long(parts[i].from, expected[2 * i]);
    test_eq_long(parts[i].to, expected[2 * i + 1]);
  }
}

void chimera_test(void) {
  Logger *logger = new_logger(stderr, LOG_ERROR);
  scheme_t *sc = scheme_new_builtin(SCHEME_DEFAULT, logger);
  enum { BODY = 200, SEG = 2 * BODY + 12 };
  const char *uni12 = "AGCAAAAGCAGG";
  const char *site = "AGCTATAGCAGG";      /* a site similar to Uni12 */
  const char *end13 = "CCTTGTTTCTACT";    /* rc of Uni13 */
  const char *start13 = "AGTAGAAACAAGG";  /* Uni13 */
  const char *gap = "ACGTTGCAACGTTGCAACGT";
  unsigned long state = 11;
  char a[BODY + 1], b[BODY + 1], other[BODY + 1], rc_a[BODY + 1];
  random_seq(a, BODY, &state);
  random_seq(b, BODY, &state);
  random_seq(other, BODY, &state);
  primers_reverse_complement(a, BODY, rc_a);
  rc_a[BODY] = '\0';
  /* the segment: two halves around the site */
  char segment[SEG + 1];
  snprintf(segment, sizeof(segment), "%s%s%s", a, site, b);
  fasta_record_t rec = {"s", "", segment, SEG};
  fasta_t fa = {&rec, 1};
  reference_index_t *ix = reference_index_new(&fa);
  reference_set_t *rs = reference_set_new(ix, 32, 0.2);
  chimera_finder_t *f = chimera_finder_new(sc, 0.2, 100);
  char read[1024];
  long whole[2] = {0, 0};

  /* one amplicon: not split */
  snprintf(read, sizeof(read), "%s%s%s", uni12, a, end13);
  whole[1] = (long)strlen(read) - 1;
  check_parts(f, read, 1, whole);

  /* two amplicons: cut after the end of the first and before the start
     of the second, on either strand */
  snprintf(read, sizeof(read), "%s%s%s%s%s%s", uni12, a, end13, gap,
           uni12, b);
  long two[4] = {0, 12 + BODY + 13 - 1, 12 + BODY + 13 + 20,
                 (long)strlen(read) - 1};
  check_parts(f, read, 2, two);
  snprintf(read, sizeof(read), "%s%s%s%s%s%s", uni12, a, end13, gap,
           start13, b);
  two[3] = (long)strlen(read) - 1;
  check_parts(f, read, 2, two);

  /* a primer inside a read, without the end of an amplicon before it:
     a junction only with references */
  snprintf(read, sizeof(read), "%s%s%s%s", uni12, other, uni12, b);
  whole[1] = (long)strlen(read) - 1;
  check_parts(f, read, 1, whole);
  chimera_finder_use_references(f, rs);
  long single[4] = {0, 12 + BODY - 1, 12 + BODY, (long)strlen(read) - 1};
  check_parts(f, read, 2, single);

  /* the end of an amplicon on the other strand, then other bases */
  snprintf(read, sizeof(read), "%s%s%s", rc_a, "CCTGCTTTTGCT", other);
  long after[4] = {0, BODY + 12 - 1, BODY + 12, (long)strlen(read) - 1};
  check_parts(f, read, 2, after);

  /* the read runs through the site of the segment: not split */
  snprintf(read, sizeof(read), "%s%s", uni12, segment);
  whole[1] = (long)strlen(read) - 1;
  check_parts(f, read, 1, whole);

  /* the same, the bases before the site differing from the reference:
     those further away still tell that the read runs through it */
  memcpy(read + 12 + BODY - 24, "ACGTACGTACGTACGTACGTACGT", 24);
  check_parts(f, read, 1, whole);

  /* an empty read */
  chimera_part_t parts[CHIMERA_MAX_PARTS];
  test_eq_size(chimera_finder_run(f, "", 0, parts), 0);

  chimera_finder_delete(f);
  reference_set_delete(rs);
  reference_index_delete(ix);
  scheme_delete(sc);
  delete_logger(logger);
}
