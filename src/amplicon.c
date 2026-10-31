/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "amplicon.h"
#include "fasta.h"
#include "logger.h"
#include "numutils.h"
#include "primer-search.h"
#include "primers.h"
#include "reference.h"
#include "scheme.h"
#include "sequtil.h"
#include "test.h"

struct amplicon_analyzer_t {
  const scheme_t *scheme;
  /* the errors allowed in a primer (see amplicon_analyzer_new()) */
  double max_errors;
  trim_policy_t policy;
  /* the primer of an end must lie within this many bases of that end of
     the read (see amplicon.h) */
  long window;
  /* the length of the longest primer of the scheme */
  size_t max_primer_len;
  /* the errors allowed in each primer of each end */
  long *left_max_errors, *right_max_errors;
  /* the costs of the missing outer bases of each primer of each end,
     max_primer_len + 1 for each primer (see missing_costs()) */
  long *left_skip, *left_edge, *right_skip, *right_edge;
  primer_search_t *search;
  /* the reverse complement of the read */
  char *rc;
  size_t rc_cap;
  /* the references, if any (see amplicon_analyzer_use_references()) */
  reference_set_t *refs;
  long terminal_window;
  /* the probe of the right end (reference_set_probe_len() + 1 bytes) */
  char *probe;
  /* the adapters, if any (see amplicon_analyzer_use_adapters()) */
  adapters_t *adapters;
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

/* the shortest truncated form of a primer searched at an end: half of
   it, at least AMPLICON_MIN_TRUNCATED bases */
static size_t min_truncated_len(size_t len) {
  size_t half = (len + 1) / 2;
  size_t min = half > AMPLICON_MIN_TRUNCATED ? half : AMPLICON_MIN_TRUNCATED;
  return min < len ? min : len;
}

/*
   The costs of the missing outer bases of a primer of length <len>, at
   an end of a read (see primer_search_end()), with the errors
   <max_errors> (see amplicon_analyzer_new()): a match is accepted if
   its errors plus the cost are at most the errors allowed in the whole
   primer (<k>), so that the cost is <k> minus the errors allowed in
   the shorter match.

   A read can begin (or end) inside a primer after a flank, e.g. when
   the amplicon was shortened before the adapter was ligated: up to
   half of the outer bases of the primer can be missing (at least
   AMPLICON_MIN_TRUNCATED bases are left), and the shorter primer,
   which matches by chance more easily, is allowed one error less than
   its length allows (<skip>). At the end of the read, the rest of the
   primer can also overhang the read, by up to half of its length, the
   errors allowed being reduced proportionally to the part inside the
   read (<edge>, the lowest cost of the two).
*/
static void missing_costs(double max_errors, size_t len, long k,
    long *skip, long *edge) {
  size_t max_t = len - min_truncated_len(len);
  for (size_t q = 0; q <= len; q++) {
    skip[q] = -1;
    edge[q] = -1;
  }
  for (size_t t = 0; t <= max_t; t++) {
    long k_t = k;
    if (t > 0) {
      k_t = max_errors_of(max_errors, len - t) - 1;
      if (k_t < 0)
        k_t = 0;
    }
    skip[t] = k - k_t;
    for (size_t o = 0; o <= (len - t) / 2; o++) {
      size_t q = t + o;
      long cost = k - floordiv(k_t * (long)(len - q), (long)(len - t));
      if (edge[q] < 0 || cost < edge[q])
        edge[q] = cost;
    }
  }
}

amplicon_analyzer_t* amplicon_analyzer_new(const scheme_t *scheme,
    double max_errors, trim_policy_t policy, long window) {
  assert(scheme != NULL);
  assert(window >= 0);
  amplicon_analyzer_t *self;
  XMALLOC(self, sizeof(amplicon_analyzer_t));
  self->scheme = scheme;
  self->max_errors = max_errors;
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
  size_t stride = self->max_primer_len + 1;
  XMALLOC(self->left_skip, sizeof(long) * stride * scheme->n_left);
  XMALLOC(self->left_edge, sizeof(long) * stride * scheme->n_left);
  XMALLOC(self->right_skip, sizeof(long) * stride * scheme->n_right);
  XMALLOC(self->right_edge, sizeof(long) * stride * scheme->n_right);
  for (size_t i = 0; i < scheme->n_left; i++)
    missing_costs(max_errors, scheme->left[i].seq_len,
                  self->left_max_errors[i], self->left_skip + i * stride,
                  self->left_edge + i * stride);
  for (size_t i = 0; i < scheme->n_right; i++)
    missing_costs(max_errors, scheme->right[i].seq_len,
                  self->right_max_errors[i], self->right_skip + i * stride,
                  self->right_edge + i * stride);
  self->search = primer_search_new();
  self->rc = NULL;
  self->rc_cap = 0;
  self->refs = NULL;
  self->terminal_window = 0;
  self->probe = NULL;
  self->adapters = NULL;
  return self;
}

void amplicon_analyzer_use_adapters(amplicon_analyzer_t *self,
    adapters_t *adapters) {
  assert(adapters != NULL);
  self->adapters = adapters_any(adapters) ? adapters : NULL;
}

void amplicon_analyzer_use_references(amplicon_analyzer_t *self,
    reference_set_t *refs, long terminal_window) {
  assert(refs != NULL);
  assert(terminal_window >= 0);
  self->refs = refs;
  self->terminal_window = terminal_window;
  XREALLOC(self->probe, reference_set_probe_len(refs) + 1);
}

void amplicon_analyzer_restrict_subtype(amplicon_analyzer_t *self,
    int subtype) {
  if (self->refs != NULL)
    reference_set_restrict(self->refs, subtype);
}

void amplicon_analyzer_delete(amplicon_analyzer_t *self) {
  if (self == NULL)
    return;
  free(self->left_max_errors);
  free(self->right_max_errors);
  free(self->left_skip);
  free(self->left_edge);
  free(self->right_skip);
  free(self->right_edge);
  primer_search_delete(self->search);
  free(self->rc);
  free(self->probe);
  free(self);
}

/*
   Whether <hit> is better than <best> as the primer of an end of a read:
   at the left end (<left>) the one ending first, at the right end the
   one starting last (the primer which started the amplicon is the
   outermost one: a primer-like sequence further inside is template);
   then the one with fewer errors, then the longer one. For a match
   inside the read (<end>: false), the one with fewer errors, then the
   longer one.
*/
static bool better_hit(const primer_hit_t *hit, const primer_hit_t *best,
    bool end, bool left) {
  if (end && hit->end != best->end && left)
    return hit->end < best->end;
  if (end && hit->start != best->start && !left)
    return hit->start > best->start;
  if (hit->n_errors != best->n_errors)
    return hit->n_errors < best->n_errors;
  return hit->end - hit->start > best->end - best->start;
}

/*
   Search the primer <p> at an end of <seq> (see primer_search_end()),
   with the costs <skip> and <edge> of its missing outer bases (see
   missing_costs()) and the errors <k> allowed in the whole primer.

   The primer may overhang the read only if it is not found otherwise
   (as a primer, whole or truncated, inside the read): the few bases of
   a primer overhanging the read match by chance more easily than the
   primer found further inside, which is then the one which started the
   amplicon.
*/
static bool search_primer_end(amplicon_analyzer_t *self, const primer_t *p,
    const char *seq, size_t len, bool left, const long *skip,
    const long *edge, long k, primer_hit_t *hit) {
  return primer_search_end(self->search, p, seq, len, left, skip, skip, k,
                           hit) ||
         primer_search_end(self->search, p, seq, len, left, skip, edge, k,
                           hit);
}

/*
   Search the primers of one end of the scheme in <seq>, and keep the
   best hit (see better_hit(); among equal hits, the first primer of the
   scheme, so that the result does not depend on the order of the
   search).

   <partial_left> and <partial_right>: whether this is the search at the
   start or at the end of the read; both are false for a match inside
   the read, which must be complete. At an end, the outer bases of a
   primer may be missing (see search_primer_end()), at the costs <skip>
   and <edge> (max_primer_len + 1 for each primer; not used for a match
   inside the read).
*/
static void search_end(amplicon_analyzer_t *self, const primer_t *primers,
    const long *max_errors, const long *skip, const long *edge, size_t n,
    const char *seq, size_t len, bool partial_left, bool partial_right,
    amplicon_end_t *out) {
  out->found = false;
  out->variant = 0;
  out->origin = ORIGIN_UNKNOWN;
  out->ref.mapped = false;
  bool end = partial_left || partial_right;
  size_t stride = self->max_primer_len + 1;
  for (size_t i = 0; i < n; i++) {
    primer_hit_t hit;
    bool found = end ?
        search_primer_end(self, primers + i, seq, len, partial_left,
                          skip + i * stride, edge + i * stride,
                          max_errors[i], &hit) :
        primer_search_run(self->search, primers + i, seq, len,
                          max_errors[i], false, false, &hit);
    if (!found)
      continue;
    if (out->found && !better_hit(&hit, &out->hit, end, partial_left))
      continue;
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
    size_t len, long left_window, long right_window, amplicon_t *out) {
  const scheme_t *sc = self->scheme;
  out->n_internal = 0;
  /* the region searched at each end: the window, plus the length of the
     longest primer, as a primer starting inside the window can end
     after it */
  long region = left_window + (long)self->max_primer_len;
  size_t left_len = (long)len < region ? len : (size_t)region;
  search_end(self, sc->left, self->left_max_errors, self->left_skip,
             self->left_edge, sc->n_left, seq, left_len, true, false,
             &out->left);
  /* the primer must begin within the window */
  if (out->left.found && out->left.hit.start > left_window)
    out->left.found = false;
  /* the right primer is searched after the left one, so that the two
     primers found belong to the same amplicon */
  long after_left = out->left.found ? out->left.hit.end + 1 : 0;
  region = right_window + (long)self->max_primer_len;
  long right_from = (long)len - region;
  if (right_from < after_left)
    right_from = after_left;
  out->right.found = false;
  out->right.variant = 0;
  out->right.origin = ORIGIN_UNKNOWN;
  out->right.ref.mapped = false;
  if (right_from < (long)len) {
    search_end(self, sc->right, self->right_max_errors, self->right_skip,
               self->right_edge, sc->n_right, seq + right_from,
               len - (size_t)right_from, false, true, &out->right);
    if (out->right.found) {
      out->right.hit.start += right_from;
      out->right.hit.end += right_from;
      out->right.hit.terminus_start += right_from;
      out->right.hit.terminus_end += right_from;
      /* the primer must end within the window */
      if (out->right.hit.end < (long)len - 1 - right_window)
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
    search_end(self, sc->left, self->left_max_errors, NULL, NULL,
               sc->n_left, seq + lo, (size_t)(hi - lo), false, false,
               &best);
    /* (the primers of the right end are searched as well: a match of
       either primer inside the read is reported) */
    amplicon_end_t other;
    search_end(self, sc->right, self->right_max_errors, NULL, NULL,
               sc->n_right, seq + lo, (size_t)(hi - lo), false, false,
               &other);
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

/* locate the probe <probe> (length <len>) of an end in the references,
   and set the origin of the end accordingly */
static void classify_end(amplicon_analyzer_t *self, const char *probe,
    size_t len, amplicon_end_t *end) {
  reference_set_locate(self->refs, probe, len, &end->ref);
  if (!end->ref.mapped)
    end->origin = ORIGIN_UNKNOWN;
  else
    end->origin = reference_hit_is_terminal(&end->ref,
        self->terminal_window) ? ORIGIN_TERMINAL : ORIGIN_INTERNAL;
}

/*
   Classify the ends of the amplicon <a> found in <seq> (length <len>)
   against the references (see amplicon.h): the probe of the left end is
   the sequence after its primer; that of the right end the reverse
   complement of the sequence before its primer, so that both are read
   from the primer into the template, as the strand the primer started.
   The probes are taken from the bases before the other primer only.
*/
static void classify_ends(amplicon_analyzer_t *self, const char *seq,
    size_t len, amplicon_t *a) {
  long lo = a->left.found ? a->left.hit.end + 1 : 0;
  long hi = a->right.found ? a->right.hit.start : (long)len;
  long avail = hi > lo ? hi - lo : 0;
  long probe_len = (long)reference_set_probe_len(self->refs);
  long n = avail < probe_len ? avail : probe_len;
  if (a->left.found)
    classify_end(self, seq + lo, (size_t)n, &a->left);
  if (a->right.found) {
    primers_reverse_complement(seq + hi - n, (size_t)n, self->probe);
    self->probe[n] = '\0';
    classify_end(self, self->probe, (size_t)n, &a->right);
  }
}

/*
   The part of <seq> (length <len>) to keep, in the orientation in which
   the layout was found; an end started by internal priming is trimmed
   as by TRIM_PRIMER, whatever the policy (see amplicon.h)
*/
static void compute_keep(const amplicon_t *a, trim_policy_t policy,
    size_t len, long *from, long *to) {
  *from = 0;
  *to = (long)len - 1;
  /* the conserved terminus is kept only when it is confidently terminal
     (TRIM_TERMINUS and ORIGIN_TERMINAL); an internal or an unlocated end
     ("cut when unsure") has its whole primer removed, so that a primer
     base never reaches the consensus unless the end is known to be a
     segment terminus */
  if (a->left.found)
    *from = policy == TRIM_TERMINUS && a->left.origin == ORIGIN_TERMINAL ?
            clamp(a->left.hit.terminus_start, 0, (long)len) :
            a->left.hit.end + 1;
  if (a->right.found)
    *to = policy == TRIM_TERMINUS && a->right.origin == ORIGIN_TERMINAL ?
          clamp(a->right.hit.terminus_end, 0, (long)len) - 1 :
          a->right.hit.start - 1;
  *from = clamp(*from, 0, (long)len);
  *to = clamp(*to, -1, (long)len - 1);
  /* no template between the two primers (a primer dimer): nothing is
     kept, whatever the policy (the termini alone are primer sequence) */
  if (a->left.found && a->right.found &&
      a->right.hit.start <= a->left.hit.end + 1) {
    *from = 1;
    *to = 0;
  }
}

/* the number of errors of the primers found (for the comparison of the
   two orientations) */
static long total_errors(const amplicon_t *a) {
  return (a->left.found ? a->left.hit.n_errors : 0) +
         (a->right.found ? a->right.hit.n_errors : 0);
}

/* the analysis of a read, without the adapters (see
   amplicon_analyzer_run()), into <out>; <start_window>, <end_window>: the
   windows at the start and at the end of the read */
static void analyze(amplicon_analyzer_t *self, const char *seq, size_t len,
    long start_window, long end_window, amplicon_t *out) {
  amplicon_t plus;
  search_layout(self, seq, len, start_window, end_window, &plus);
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
    search_layout(self, self->rc, len, end_window, start_window, &minus);
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
  if (self->refs != NULL)
    classify_ends(self, use_minus ? self->rc : seq, len, out);
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
}

/* move the positions of <a>, the analysis of the part of a read of
   length <len> from <from> to <to>, to the whole read */
static void shift_positions(amplicon_t *a, long from, long to, size_t len) {
  if (a->keep_to >= a->keep_from) {
    a->keep_from += from;
    a->keep_to += from;
  }
  for (size_t i = 0; i < a->n_internal; i++) {
    a->internal[i].start += from;
    a->internal[i].end += from;
  }
  /* the primers: in the analyzed orientation (for the minus strand, the
     reverse complement of the read) */
  long shift = a->strand == STRAND_MINUS ? (long)len - 1 - to : from;
  amplicon_end_t *ends[] = {&a->left, &a->right};
  for (size_t e = 0; e < 2; e++)
    if (ends[e]->found) {
      ends[e]->hit.start += shift;
      ends[e]->hit.end += shift;
      ends[e]->hit.terminus_start += shift;
      ends[e]->hit.terminus_end += shift;
    }
}

const amplicon_t* amplicon_analyzer_run(amplicon_analyzer_t *self,
    const char *seq, size_t len) {
  assert(self != NULL);
  assert(seq != NULL);
  amplicon_t *out = &self->result;
  analyze(self, seq, len, self->window, self->window, out);
  out->adapter_start.found = out->adapter_end.found = false;
  if (self->adapters == NULL)
    return out;
  /* the adapters, at the ends of the read without a primer */
  bool primer_start = out->strand == STRAND_MINUS ? out->right.found :
                      out->left.found;
  bool primer_end = out->strand == STRAND_MINUS ? out->left.found :
                    out->right.found;
  adapter_hit_t a_start = {.found = false}, a_end = {.found = false};
  if (!primer_start)
    adapters_search_start(self->adapters, seq, len, &a_start);
  if (!primer_end)
    adapters_search_end(self->adapters, seq, len, &a_end);
  if (!a_start.found && !a_end.found)
    return out;
  long from = a_start.found ? a_start.end + 1 : 0;
  long to = a_end.found ? a_end.start - 1 : (long)len - 1;
  if (to < from) {
    /* nothing between the adapters */
    analyze(self, seq, 0, self->window, self->window, out);
    out->keep_from = 1;
    out->keep_to = 0;
  } else {
    /* the part between the adapters is analyzed again: a primer whose
       outer bases are missing is now at an end, next to the adapter (see
       amplicon.h) */
    amplicon_t whole = *out;
    analyze(self, seq + from, (size_t)(to - from + 1),
            a_start.found ? AMPLICON_ADAPTER_WINDOW : self->window,
            a_end.found ? AMPLICON_ADAPTER_WINDOW : self->window, out);
    if (amplicon_n_ends(out) < amplicon_n_ends(&whole))
      *out = whole;
    else
      shift_positions(out, from, to, len);
    if (out->keep_from < from)
      out->keep_from = from;
    if (out->keep_to > to)
      out->keep_to = to;
  }
  out->adapter_start = a_start;
  out->adapter_end = a_end;
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

/* the reverse complement of <seq> into <dest> (0-terminated) */
static void test_rc(const char *seq, char *dest) {
  size_t len = strlen(seq);
  reverse_complement(seq, len, dest);
  dest[len] = '\0';
}

/*
   The classification of the ends against the references, with the
   default scheme and policy: the conserved termini of the reads of a
   normal amplicon are kept, the primer of an internally primed read is
   removed whole.
*/
static void test_references(Logger *logger) {
  /* the interior of a segment (not similar to the primers) */
  const char *interior =
    "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG"
    "CTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGAC"
    "TGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGACAGGTCACGCAGAGGC"
    "GCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAA"
    "CTCCAGCGCGGTCAGTTCCATCACCCTAAGTAACCGAATAATGCGTTCGCTCTATTGACT";
  const char *uni12 = "AGCAAAAGCAGG";
  const char *end13 = "CCTTGTTTCTACT";                /* rc of Uni13 */
  const char *rc13 = "CCTTGTTTCTACTGATCACGCGT";       /* rc mbtuni-13 */
  /* the segment: the conserved termini and the interior */
  char segment[512];
  snprintf(segment, sizeof(segment), "%s%s%s", uni12, interior, end13);
  fasta_record_t record = {.name = NULL, .description = NULL,
                           .seq = segment, .seq_len = strlen(segment)};
  fasta_t fa = {.records = &record, .n_records = 1};
  reference_index_t *index = reference_index_new(&fa);
  reference_set_t *refs = reference_set_new(index,
      REFERENCE_DEFAULT_PROBE_LEN, REFERENCE_DEFAULT_MAX_ERROR_RATE);
  scheme_t *sc = scheme_new_builtin(SCHEME_DEFAULT, logger);
  amplicon_analyzer_t *a = amplicon_analyzer_new(sc, 0.1, TRIM_TERMINUS,
      AMPLICON_DEFAULT_WINDOW);
  amplicon_analyzer_use_references(a, refs,
                                   REFERENCE_DEFAULT_TERMINAL_WINDOW);
  char read[512], rc[512], expected[512], rc_expected[512];

  /* a normal amplicon: both ends terminal, the termini are kept and
     only the tail of the Uni13 primer is removed */
  snprintf(read, sizeof(read), "%s%s%s", uni12, interior, rc13);
  const amplicon_t *r = amplicon_analyzer_run(a, read, strlen(read));
  test_assert(r->left.origin == ORIGIN_TERMINAL);
  test_assert(r->right.origin == ORIGIN_TERMINAL);
  test_eq_long(r->left.ref.ref_start, (long)strlen(uni12));
  test_eq_mem(read + r->keep_from, amplicon_kept_length(r), segment);

  /* internal priming: the Uni12 primer started the amplicon 150 bases
     into the segment; its bases are not the template there and are
     removed, the rest of the read is kept */
  const char *site = interior + 150;
  snprintf(read, sizeof(read), "%s%s%s", uni12, site, rc13);
  r = amplicon_analyzer_run(a, read, strlen(read));
  test_assert(r->strand == STRAND_PLUS);
  test_assert(r->left.origin == ORIGIN_INTERNAL);
  test_assert(r->right.origin == ORIGIN_TERMINAL);
  test_eq_long(r->left.ref.ref_start, (long)strlen(uni12) + 150);
  snprintf(expected, sizeof(expected), "%s%s", site, end13);
  test_eq_mem(read + r->keep_from, amplicon_kept_length(r), expected);

  /* the same read sequenced in the other orientation: the primer is at
     the end of the original read, and is removed there */
  test_rc(read, rc);
  test_rc(expected, rc_expected);
  r = amplicon_analyzer_run(a, rc, strlen(rc));
  test_assert(r->strand == STRAND_MINUS);
  test_assert(r->left.origin == ORIGIN_INTERNAL);
  test_eq_mem(rc + r->keep_from, amplicon_kept_length(r), rc_expected);

  /* a read starting at an interior site without the right primer (a
     partial read) is classified all the same */
  snprintf(read, sizeof(read), "%s%.100s", uni12, site);
  r = amplicon_analyzer_run(a, read, strlen(read));
  test_eq_size(amplicon_n_ends(r), 1);
  test_assert(r->left.origin == ORIGIN_INTERNAL);
  snprintf(expected, sizeof(expected), "%.100s", site);
  test_eq_mem(read + r->keep_from, amplicon_kept_length(r), expected);

  /* a template absent from the references: not classified (unknown), so
     cut when unsure (the whole primer removed) */
  snprintf(read, sizeof(read), "%s%s", uni12,
           "GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAA");
  r = amplicon_analyzer_run(a, read, strlen(read));
  test_assert(r->left.found);
  test_assert(r->left.origin == ORIGIN_UNKNOWN);
  test_assert(!r->left.ref.mapped);
  test_eq_long(r->keep_from, (long)strlen(uni12));

  /* a primer dimer: no probe to locate */
  snprintf(read, sizeof(read), "%s%s", uni12, rc13);
  r = amplicon_analyzer_run(a, read, strlen(read));
  test_assert(r->left.origin == ORIGIN_UNKNOWN);
  test_assert(r->right.origin == ORIGIN_UNKNOWN);
  amplicon_analyzer_delete(a);

  /* with the default policy, every primer is removed: the internally
     primed read loses its primer with or without the classification,
     which is still reported */
  a = amplicon_analyzer_new(sc, 0.1, AMPLICON_DEFAULT_POLICY,
                            AMPLICON_DEFAULT_WINDOW);
  amplicon_analyzer_use_references(a, refs,
                                   REFERENCE_DEFAULT_TERMINAL_WINDOW);
  snprintf(read, sizeof(read), "%s%s%s", uni12, site, rc13);
  r = amplicon_analyzer_run(a, read, strlen(read));
  test_assert(r->left.origin == ORIGIN_INTERNAL);
  test_eq_mem(read + r->keep_from, amplicon_kept_length(r), site);
  snprintf(read, sizeof(read), "%s%s%s", uni12, interior, rc13);
  r = amplicon_analyzer_run(a, read, strlen(read));
  test_assert(r->left.origin == ORIGIN_TERMINAL);
  test_eq_mem(read + r->keep_from, amplicon_kept_length(r), interior);
  amplicon_analyzer_delete(a);

  /* TRIM_TERMINUS without the references: the internally primed end is
     unknown (no classification) and cut when unsure, so its primer is
     removed and the artefact cannot reach the consensus */
  a = amplicon_analyzer_new(sc, 0.1, TRIM_TERMINUS,
                            AMPLICON_DEFAULT_WINDOW);
  snprintf(read, sizeof(read), "%s%s%s", uni12, site, rc13);
  r = amplicon_analyzer_run(a, read, strlen(read));
  test_assert(r->left.origin == ORIGIN_UNKNOWN);
  test_eq_long(r->keep_from, (long)strlen(uni12));
  amplicon_analyzer_delete(a);
  scheme_delete(sc);
  reference_set_delete(refs);
  reference_index_delete(index);
}

/* the adapters at the ends without a primer (see amplicon.h) */
static void test_adapters(const scheme_t *sc, const char *body,
    const char *p12, const char *rc13) {
  /* native barcode NB01: the adapter at the start of a read, and the
     reverse complement of the bottom adapter at its end */
  const char *nb_start = "ATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCT";
  const char *nb_end = "GGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCTTAGCAAT";
  amplicon_analyzer_t *a = amplicon_analyzer_new(sc, 0.2, TRIM_PRIMER,
      AMPLICON_DEFAULT_WINDOW);
  adapters_t *ad = adapters_new("nb", AMPLICON_DEFAULT_WINDOW, NULL);
  amplicon_analyzer_use_adapters(a, ad);
  char read[512], rc[512];
  /* a broken read: no primer at the start, the adapter is removed */
  snprintf(read, sizeof(read), "GGCAT%s%s%s", nb_start, body, rc13);
  check_keep(a, read, STRAND_PLUS, 1, body);
  const amplicon_t *r = amplicon_analyzer_run(a, read, strlen(read));
  test_assert(r->adapter_start.found);
  test_assert(!r->adapter_end.found);
  /* a broken read of the other strand (the adapters at the same ends
     of the read): no primer at the end */
  char amp[256];
  snprintf(amp, sizeof(amp), "%s%s", body, rc13);
  test_rc(amp, rc);
  snprintf(read, sizeof(read), "%s%s%sGA", nb_start, rc, nb_end);
  r = amplicon_analyzer_run(a, read, strlen(read));
  test_assert(r->strand == STRAND_MINUS);
  test_assert(!r->adapter_start.found);
  test_assert(r->adapter_end.found);
  test_eq_size(amplicon_kept_length(r), strlen(body));
  test_eq_long(r->keep_from, (long)(strlen(nb_start) + strlen(rc13)));
  /* a primer without its first 14 bases after the adapter: too short to
     be searched after a flank, removed next to the adapter */
  snprintf(read, sizeof(read), "%s%s%s%s", nb_start, p12 + 14, body, rc13);
  check_keep(a, read, STRAND_PLUS, 2, body);
  /* an adapter at the end, with no primer before it */
  snprintf(read, sizeof(read), "%s%s%sTTAG", p12, body, nb_end);
  check_keep(a, read, STRAND_PLUS, 1, body);
  /* no primer and no adapter: nothing is removed */
  check_keep(a, body, STRAND_NONE, 0, body);
  /* the primers found: the adapters are not searched */
  snprintf(read, sizeof(read), "%s%s%s%s%s", nb_start, p12, body, rc13,
           nb_end);
  check_keep(a, read, STRAND_PLUS, 2, body);
  r = amplicon_analyzer_run(a, read, strlen(read));
  test_assert(!r->adapter_start.found && !r->adapter_end.found);
  /* only adapters: nothing is kept */
  snprintf(read, sizeof(read), "%s%s", nb_start, nb_end);
  r = amplicon_analyzer_run(a, read, strlen(read));
  test_eq_size(amplicon_kept_length(r), 0);
  adapters_delete(ad);
  amplicon_analyzer_delete(a);
}

void amplicon_test(void) {
  Logger *logger = new_logger(stderr, LOG_ERROR);
  /* the costs of the missing outer bases of a primer of 12 bases with
     errors 0.2 (two errors in the whole primer): one or two bases
     missing after a flank leave one error (cost 1), three or four
     none (cost 2); at the end of the read, up to six bases can overhang
     with one error, and the truncated primers can overhang as well */
  long skip[13], edge[13];
  const long exp_skip[13] = {0, 1, 1, 2, 2, -1, -1, -1, -1, -1, -1, -1, -1};
  const long exp_edge[13] = {0, 1, 1, 1, 1, 1, 1, 2, 2, -1, -1, -1, -1};
  missing_costs(0.2, 12, 2, skip, edge);
  for (size_t q = 0; q <= 12; q++) {
    test_eq_long(skip[q], exp_skip[q]);
    test_eq_long(edge[q], exp_edge[q]);
  }
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

  /* TRIM_TERMINUS without references: no end can be classified as
     terminal, so every end is cut when unsure (the whole primer
     removed, as in TRIM_PRIMER); keeping the termini needs the
     references (see test_references) */
  a = amplicon_analyzer_new(sc, 0.2, TRIM_TERMINUS,
                            AMPLICON_DEFAULT_WINDOW);
  snprintf(read, sizeof(read), "%s%s%s", p12, body, rc13);
  check_keep(a, read, STRAND_PLUS, 2, body);
  /* a primer dimer: nothing is kept, also with this policy */
  snprintf(read, sizeof(read), "%s%s", p12, rc13);
  r = amplicon_analyzer_run(a, read, strlen(read));
  test_eq_size(amplicon_n_ends(r), 2);
  test_eq_size(amplicon_kept_length(r), 0);
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

  /*
     A read which begins with flanking bases and then inside the primer
     (here without its first two bases), and has a primer-like sequence
     inside the window (the site of the HA of KF918483.1 at 1145-1157,
     edit distance 2 from Uni12): the truncated primer, outermost, is
     the primer of the left end, and the template is kept from after it
  */
  a = amplicon_analyzer_new(sc, 0.2, TRIM_PRIMER, AMPLICON_DEFAULT_WINDOW);
  char trunc[512];
  snprintf(trunc, sizeof(trunc), "TTGCA%s%sAGCAACGAGCAGG%s%s",
           uni12 + 2, body, body, rc13);
  r = amplicon_analyzer_run(a, trunc, strlen(trunc));
  test_eq(r->keep_from, (long)(5 + strlen(uni12) - 2), long, "%ld");
  /* without the truncated primer, the primer-like sequence would be
     taken: it is the only match of the whole primer */
  snprintf(trunc, sizeof(trunc), "TTGCAGG%sAGCAACGAGCAGG%s%s", body,
           body, rc13);
  r = amplicon_analyzer_run(a, trunc, strlen(trunc));
  test_eq(r->keep_from, (long)(7 + strlen(body) + 13), long, "%ld");
  /* the truncated primer is allowed one error less than its length
     allows: with one mismatch, it is still taken (its cost equals that
     of the primer-like sequence, and it is the outermost); with two, the
     primer-like sequence is taken */
  snprintf(trunc, sizeof(trunc), "TTGCACAAATGCAGG%sAGCAACGAGCAGG%s%s",
           body, body, rc13);
  r = amplicon_analyzer_run(a, trunc, strlen(trunc));
  test_eq(r->keep_from, 15, long, "%ld");
  snprintf(trunc, sizeof(trunc), "TTGCACAAATGCTGG%sAGCAACGAGCAGG%s%s",
           body, body, rc13);
  r = amplicon_analyzer_run(a, trunc, strlen(trunc));
  test_eq(r->keep_from, (long)(15 + strlen(body) + 13), long, "%ld");
  /* the primer of the right end without its last five bases, followed
     by flanking bases */
  snprintf(trunc, sizeof(trunc), "%s%s%.18sAAGTTC", uni12, body, rc13);
  check_keep(a, trunc, STRAND_PLUS, 2, body);
  amplicon_analyzer_delete(a);

  /* with a window wide enough, the same match is taken as the primer of
     the left end (which is why the window matters) */
  a = amplicon_analyzer_new(sc, 0.1, TRIM_PRIMER, 1000);
  r = amplicon_analyzer_run(a, internal_only, strlen(internal_only));
  test_eq_size(amplicon_n_ends(r), 1);
  test_eq(r->keep_from, (long)(4 * strlen(body) + strlen(uni12)), long,
          "%ld");
  amplicon_analyzer_delete(a);
  test_adapters(sc, body, p12, rc13);
  scheme_delete(sc);
  test_references(logger);
  delete_logger(logger);
}
