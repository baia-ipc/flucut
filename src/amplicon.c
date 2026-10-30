/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "amplicon.h"
#include "logger.h"
#include "primer-search.h"
#include "primers.h"
#include "scheme.h"
#include "sequtil.h"
#include "test.h"

struct amplicon_analyzer_t {
  const scheme_t *scheme;
  trim_policy_t policy;
  /* the primer of an end must lie within this many bases of that end of
     the read (see amplicon.h) */
  long window;
  /* the length of the longest primer of the scheme */
  size_t max_primer_len;
  /* the errors allowed in each primer of each end */
  long *left_max_errors, *right_max_errors;
  primer_search_t *search;
  /* the reverse complement of the read */
  char *rc;
  size_t rc_cap;
  amplicon_t result;
};

/* the errors allowed in a primer: a number of bases, or a fraction of
   its length if smaller than 1 */
static long max_errors_of(double max_errors, size_t seq_len) {
  if (max_errors >= 1.0)
    return (long)max_errors;
  long n = (long)(max_errors * (double)seq_len);
  return n > 0 ? n : 0;
}

amplicon_analyzer_t* amplicon_analyzer_new(const scheme_t *scheme,
    double max_errors, trim_policy_t policy, long window) {
  assert(scheme != NULL);
  assert(window >= 0);
  amplicon_analyzer_t *self;
  XMALLOC(self, sizeof(amplicon_analyzer_t));
  self->scheme = scheme;
  self->policy = policy;
  self->window = window;
  self->max_primer_len = 0;
  for (size_t i = 0; i < scheme->n_left; i++)
    if (scheme->left[i].seq_len > self->max_primer_len)
      self->max_primer_len = scheme->left[i].seq_len;
  for (size_t i = 0; i < scheme->n_right; i++)
    if (scheme->right[i].seq_len > self->max_primer_len)
      self->max_primer_len = scheme->right[i].seq_len;
  XMALLOC(self->left_max_errors, sizeof(long) * scheme->n_left);
  XMALLOC(self->right_max_errors, sizeof(long) * scheme->n_right);
  for (size_t i = 0; i < scheme->n_left; i++)
    self->left_max_errors[i] = max_errors_of(max_errors,
                                             scheme->left[i].seq_len);
  for (size_t i = 0; i < scheme->n_right; i++)
    self->right_max_errors[i] = max_errors_of(max_errors,
                                              scheme->right[i].seq_len);
  self->search = primer_search_new();
  self->rc = NULL;
  self->rc_cap = 0;
  return self;
}

void amplicon_analyzer_delete(amplicon_analyzer_t *self) {
  if (self == NULL)
    return;
  free(self->left_max_errors);
  free(self->right_max_errors);
  primer_search_delete(self->search);
  free(self->rc);
  free(self);
}

/*
   Search the primers of one end of the scheme in <seq>, and keep the
   best hit: the one with the fewest errors (among those, the longest
   alignment, then the first primer of the scheme, so that the result
   does not depend on the order of the search).

   <partial_left> and <partial_right>: whether the primer may overhang
   that side of <seq>; both are false for a match inside the read, which
   must be complete.
*/
static void search_end(amplicon_analyzer_t *self, const primer_t *primers,
    const long *max_errors, size_t n, const char *seq, size_t len,
    bool partial_left, bool partial_right, amplicon_end_t *out) {
  out->found = false;
  out->variant = 0;
  for (size_t i = 0; i < n; i++) {
    primer_hit_t hit;
    if (!primer_search_run(self->search, primers + i, seq, len,
                           max_errors[i], partial_left, partial_right,
                           &hit))
      continue;
    if (out->found) {
      long best = out->hit.end - out->hit.start;
      long cur = hit.end - hit.start;
      if (hit.n_errors > out->hit.n_errors ||
          (hit.n_errors == out->hit.n_errors && cur <= best))
        continue;
    }
    out->found = true;
    out->variant = i;
    out->hit = hit;
  }
}

/*
   Search the layout of the plus strand of an amplicon in <seq>: the
   primers of the left end at its start, those of the right end at its
   end.

   Each primer is searched only in the window at its end of the read
   (see amplicon.h): a match further inside the read is not the primer
   which started the amplicon, and trimming at it would remove template
   sequence. Such matches are collected separately, as internal matches.
*/
static void search_layout(amplicon_analyzer_t *self, const char *seq,
    size_t len, amplicon_t *out) {
  const scheme_t *sc = self->scheme;
  out->n_internal = 0;
  /* the region searched at each end: the window, plus the length of the
     longest primer, as a primer starting inside the window can end
     after it */
  long region = self->window + (long)self->max_primer_len;
  size_t left_len = (long)len < region ? len : (size_t)region;
  search_end(self, sc->left, self->left_max_errors, sc->n_left, seq,
             left_len, true, false, &out->left);
  /* the primer must begin within the window */
  if (out->left.found && out->left.hit.start > self->window)
    out->left.found = false;
  /* the right primer is searched after the left one, so that the two
     primers found belong to the same amplicon */
  long after_left = out->left.found ? out->left.hit.end + 1 : 0;
  long right_from = (long)len - region;
  if (right_from < after_left)
    right_from = after_left;
  out->right.found = false;
  out->right.variant = 0;
  if (right_from < (long)len) {
    search_end(self, sc->right, self->right_max_errors, sc->n_right,
               seq + right_from, len - (size_t)right_from, false, true,
               &out->right);
    if (out->right.found) {
      out->right.hit.start += right_from;
      out->right.hit.end += right_from;
      out->right.hit.terminus_start += right_from;
      out->right.hit.terminus_end += right_from;
      /* the primer must end within the window */
      if (out->right.hit.end < (long)len - 1 - self->window)
        out->right.found = false;
    }
  }
}

/*
   Collect the matches of the primers of the scheme which lie between the
   primers of the two ends, i.e. inside the part of the read which is
   kept (see amplicon.h: internal priming and chimeras).

   The matches are searched from left to right, without a partial
   alignment at either side, and are not trimmed.
*/
static void search_internal(amplicon_analyzer_t *self, const char *seq,
    size_t len, amplicon_t *out) {
  const scheme_t *sc = self->scheme;
  long lo = out->left.found ? out->left.hit.end + 1 : 0;
  long hi = out->right.found ? out->right.hit.start : (long)len;
  while (out->n_internal < AMPLICON_MAX_INTERNAL && hi - lo > 0) {
    amplicon_end_t best;
    bool from_left_end = true;
    search_end(self, sc->left, self->left_max_errors, sc->n_left,
               seq + lo, (size_t)(hi - lo), false, false, &best);
    /* (the primers of the right end are searched as well: a match of
       either primer inside the read is reported) */
    amplicon_end_t other;
    search_end(self, sc->right, self->right_max_errors, sc->n_right,
               seq + lo, (size_t)(hi - lo), false, false, &other);
    if (other.found && (!best.found ||
                        other.hit.n_errors < best.hit.n_errors)) {
      best = other;
      from_left_end = false;
    }
    if (!best.found)
      return;
    amplicon_match_t *m = out->internal + out->n_internal++;
    m->at_left_end = from_left_end;
    m->variant = best.variant;
    m->start = best.hit.start + lo;
    m->end = best.hit.end + lo;
    m->n_errors = best.hit.n_errors;
    lo = m->end + 1;
  }
}

/* a value limited to [lo, hi] */
static long clamp(long v, long lo, long hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}

/*
   The part of <seq> (length <len>) to keep, in the orientation in which
   the layout was found
*/
static void compute_keep(const amplicon_t *a, trim_policy_t policy,
    size_t len, long *from, long *to) {
  *from = 0;
  *to = (long)len - 1;
  if (a->left.found)
    *from = policy == TRIM_PRIMER ? a->left.hit.end + 1 :
            clamp(a->left.hit.terminus_start, 0, (long)len);
  if (a->right.found)
    *to = policy == TRIM_PRIMER ? a->right.hit.start - 1 :
          clamp(a->right.hit.terminus_end, 0, (long)len) - 1;
  *from = clamp(*from, 0, (long)len);
  *to = clamp(*to, -1, (long)len - 1);
}

/* the number of errors of the primers found (for the comparison of the
   two orientations) */
static long total_errors(const amplicon_t *a) {
  return (a->left.found ? a->left.hit.n_errors : 0) +
         (a->right.found ? a->right.hit.n_errors : 0);
}

const amplicon_t* amplicon_analyzer_run(amplicon_analyzer_t *self,
    const char *seq, size_t len) {
  assert(self != NULL);
  assert(seq != NULL);
  amplicon_t *out = &self->result;
  amplicon_t plus;
  search_layout(self, seq, len, &plus);
  amplicon_t minus;
  bool minus_searched = false;
  /*
     The primers of the two ends of a scheme can be similar to each
     other: the tailed primers share their tail, and only their
     conserved termini differ. A primer of one end can therefore be
     found, with errors, at the other end of a read of the opposite
     orientation. Both orientations are thus searched and compared,
     except when the plus strand already gives a complete layout
     without any error, which no other orientation can better.
  */
  if (amplicon_n_ends(&plus) < 2 || total_errors(&plus) > 0) {
    XCHECKSPACEINC(self->rc, &self->rc_cap, len + 1, 1024);
    reverse_complement(seq, len, self->rc);
    self->rc[len] = '\0';
    search_layout(self, self->rc, len, &minus);
    minus_searched = true;
  }
  /* the orientation with more primers found; with the same number, the
     one with fewer errors; with the same number of errors, the plus
     strand */
  bool use_minus = minus_searched && amplicon_n_ends(&minus) > 0 &&
      (amplicon_n_ends(&minus) > amplicon_n_ends(&plus) ||
       (amplicon_n_ends(&minus) == amplicon_n_ends(&plus) &&
        total_errors(&minus) < total_errors(&plus)));
  *out = use_minus ? minus : plus;
  search_internal(self, use_minus ? self->rc : seq, len, out);
  out->strand = amplicon_n_ends(out) == 0 ? STRAND_NONE :
                (use_minus ? STRAND_MINUS : STRAND_PLUS);
  long from, to;
  compute_keep(out, self->policy, len, &from, &to);
  /* the positions of the original read: for a read of the minus strand,
     the kept part is mirrored back (the read is not reverse-complemented,
     see amplicon.h) */
  if (use_minus) {
    long n = (long)len;
    out->keep_from = to < from ? 1 : n - 1 - to;
    out->keep_to = to < from ? 0 : n - 1 - from;
    for (size_t i = 0; i < out->n_internal; i++) {
      amplicon_match_t *m = out->internal + i;
      long start = n - 1 - m->end;
      m->end = n - 1 - m->start;
      m->start = start;
    }
  } else {
    out->keep_from = from;
    out->keep_to = to;
  }
  return out;
}

/* analyze a read and check the kept part */
static void check_keep(amplicon_analyzer_t *a, const char *read,
    strand_t strand, size_t n_ends, const char *kept) {
  const amplicon_t *r = amplicon_analyzer_run(a, read, strlen(read));
  test_assert(r->strand == strand);
  test_eq_size(amplicon_n_ends(r), n_ends);
  test_eq_mem(read + r->keep_from, amplicon_kept_length(r), kept);
}

void amplicon_test(void) {
  Logger *logger = new_logger(stderr, LOG_ERROR);
  /* the scheme with a tail at both ends */
  scheme_t *sc = scheme_new_builtin("mbt", logger);
  /* a body of 40 bases, which is not similar to the primers */
  const char *body = "GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAA";
  const char *p12 = "ACGCGTGATCAGCAAAAGCAGG";          /* mbtuni-12 */
  const char *rc13 = "CCTTGTTTCTACTGATCACGCGT";        /* rc mbtuni-13 */
  char read[256], rc[256];
  snprintf(read, sizeof(read), "%s%s%s", p12, body, rc13);

  /* the default: all the primer-derived bases are removed */
  amplicon_analyzer_t *a = amplicon_analyzer_new(sc, 0.2, TRIM_PRIMER,
      AMPLICON_DEFAULT_WINDOW);
  check_keep(a, read, STRAND_PLUS, 2, body);
  /* the same read, sequenced in the other orientation: the kept part is
     the reverse complement of the body, in the original orientation */
  reverse_complement(read, strlen(read), rc);
  rc[strlen(read)] = '\0';
  char rc_body[64];
  reverse_complement(body, strlen(body), rc_body);
  rc_body[strlen(body)] = '\0';
  check_keep(a, rc, STRAND_MINUS, 2, rc_body);
  /* a read which begins inside the left primer: what is left of the
     primer is still removed */
  snprintf(read, sizeof(read), "%s%s%s", p12 + 11, body, rc13);
  check_keep(a, read, STRAND_PLUS, 2, body);
  /* a read without the right primer: only the left end is trimmed */
  snprintf(read, sizeof(read), "%s%s", p12, body);
  check_keep(a, read, STRAND_PLUS, 1, body);
  /* a read with no primer at all is not trimmed */
  check_keep(a, body, STRAND_NONE, 0, body);
  /* a primer dimer: nothing is left */
  snprintf(read, sizeof(read), "%s%s", p12, rc13);
  const amplicon_t *r = amplicon_analyzer_run(a, read, strlen(read));
  test_eq_size(amplicon_n_ends(r), 2);
  test_eq_size(amplicon_kept_length(r), 0);
  amplicon_analyzer_delete(a);
  r = NULL;

  /* keeping the conserved termini: the body plus the 12 and 13 bases of
     the termini, which come from the primers */
  a = amplicon_analyzer_new(sc, 0.2, TRIM_TERMINUS,
                            AMPLICON_DEFAULT_WINDOW);
  snprintf(read, sizeof(read), "%s%s%s", p12, body, rc13);
  char expected[256];
  snprintf(expected, sizeof(expected), "%s%s%s", PRIMERS_TERMINUS_12,
           body, "CCTTGTTTCTACT");
  check_keep(a, read, STRAND_PLUS, 2, expected);
  amplicon_analyzer_delete(a);
  scheme_delete(sc);

  /*
     Internal priming, with the default scheme: the untailed Uni12
     primer is only the twelve bases of the conserved terminus and can
     anneal inside a segment.
  */
  sc = scheme_new_builtin(SCHEME_DEFAULT, logger);
  a = amplicon_analyzer_new(sc, 0.1, TRIM_PRIMER,
                            AMPLICON_DEFAULT_WINDOW);
  const char *uni12 = "AGCAAAAGCAGG";
  /*
     A read started at an internal site looks like a read of a real
     amplicon: the primer is at its end, so its primer-derived bases are
     removed as well and never reach the consensus
  */
  snprintf(read, sizeof(read), "%s%s%s", uni12, body, rc13);
  check_keep(a, read, STRAND_PLUS, 2, body);

  /*
     A primer which matches deep inside a read is not the primer which
     started the amplicon: it is reported, and nothing is trimmed there
  */
  char middle[256];
  snprintf(middle, sizeof(middle), "%s%s%s%s%s", uni12, body, uni12, body,
           rc13);
  r = amplicon_analyzer_run(a, middle, strlen(middle));
  test_eq_size(amplicon_n_ends(r), 2);
  test_eq_size(r->n_internal, 1);
  test_assert(r->internal[0].at_left_end);
  test_eq(r->internal[0].start, (long)(strlen(uni12) + strlen(body)),
          long, "%ld");
  /* the template between the two end primers is kept whole, including
     the internal match */
  test_eq_size(amplicon_kept_length(r),
               strlen(body) * 2 + strlen(uni12));

  /*
     Without the window, a read whose only primer match is deep inside
     it would be trimmed there, removing real template. With the
     window, no end primer is found and the read is not trimmed.
  */
  char internal_only[512];
  snprintf(internal_only, sizeof(internal_only),
           "%s%s%s%s%s%s%s%s%s", body, body, body, body, uni12, body,
           body, body, body);
  check_keep(a, internal_only, STRAND_NONE, 0, internal_only);
  amplicon_analyzer_delete(a);

  /* with a window wide enough, the same match is taken as the primer of
     the left end (which is why the window matters) */
  a = amplicon_analyzer_new(sc, 0.1, TRIM_PRIMER, 1000);
  r = amplicon_analyzer_run(a, internal_only, strlen(internal_only));
  test_eq_size(amplicon_n_ends(r), 1);
  test_eq(r->keep_from, (long)(4 * strlen(body) + strlen(uni12)), long,
          "%ld");
  amplicon_analyzer_delete(a);
  scheme_delete(sc);
  delete_logger(logger);
}
