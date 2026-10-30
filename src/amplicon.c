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
    double max_errors, trim_policy_t policy) {
  assert(scheme != NULL);
  amplicon_analyzer_t *self;
  XMALLOC(self, sizeof(amplicon_analyzer_t));
  self->scheme = scheme;
  self->policy = policy;
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

   <at_left>: the end searched is the left one, so the primer may
   overhang the start of the read; otherwise it may overhang its end.
*/
static void search_end(amplicon_analyzer_t *self, const primer_t *primers,
    const long *max_errors, size_t n, const char *seq, size_t len,
    bool at_left, amplicon_end_t *out) {
  out->found = false;
  out->variant = 0;
  for (size_t i = 0; i < n; i++) {
    primer_hit_t hit;
    if (!primer_search_run(self->search, primers + i, seq, len,
                           max_errors[i], at_left, !at_left, &hit))
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
   end. The right end is searched only after the left one, in the part
   of the read which follows it, so that the two primers found belong to
   the same amplicon.
*/
static void search_layout(amplicon_analyzer_t *self, const char *seq,
    size_t len, amplicon_t *out) {
  const scheme_t *sc = self->scheme;
  search_end(self, sc->left, self->left_max_errors, sc->n_left, seq, len,
             true, &out->left);
  /* the right primer is searched after the left one: with a left primer
     found, the search starts after its last base */
  long from = out->left.found ? out->left.hit.end + 1 : 0;
  if (from >= (long)len) {
    out->right.found = false;
    out->right.variant = 0;
    return;
  }
  search_end(self, sc->right, self->right_max_errors, sc->n_right,
             seq + from, len - (size_t)from, false, &out->right);
  if (out->right.found) {
    out->right.hit.start += from;
    out->right.hit.end += from;
    out->right.hit.terminus_start += from;
    out->right.hit.terminus_end += from;
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
  scheme_t *sc = scheme_new_builtin(SCHEME_DEFAULT, logger);
  /* a body of 40 bases, which is not similar to the primers */
  const char *body = "GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAA";
  const char *p12 = "ACGCGTGATCAGCAAAAGCAGG";          /* mbtuni-12 */
  const char *rc13 = "CCTTGTTTCTACTGATCACGCGT";        /* rc mbtuni-13 */
  char read[256], rc[256];
  snprintf(read, sizeof(read), "%s%s%s", p12, body, rc13);

  /* the default: all the primer-derived bases are removed */
  amplicon_analyzer_t *a = amplicon_analyzer_new(sc, 0.2, TRIM_PRIMER);
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

  /* keeping the conserved termini: the body plus the 12 and 13 bases of
     the termini, which come from the primers */
  a = amplicon_analyzer_new(sc, 0.2, TRIM_TERMINUS);
  snprintf(read, sizeof(read), "%s%s%s", p12, body, rc13);
  char expected[256];
  snprintf(expected, sizeof(expected), "%s%s%s", PRIMERS_TERMINUS_12,
           body, "CCTTGTTTCTACT");
  check_keep(a, read, STRAND_PLUS, 2, expected);
  amplicon_analyzer_delete(a);
  scheme_delete(sc);
  delete_logger(logger);
}
