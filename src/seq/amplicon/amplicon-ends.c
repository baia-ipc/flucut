/* (c) 2026 Giorgio Gonnella */

#include <stdbool.h>
#include <stddef.h>
#include "amplicon-analyzer.h"
#include "amplicon-ends.h"
#include "internal-matches.h"
#include "numutils.h"
#include "primer-search.h"
#include "primers.h"
#include "scheme.h"
#include "virus.h"

/* the errors allowed in a primer: a number of bases, or a fraction of
   its length if smaller than 1 */
long amplicon_max_errors_of(double max_errors, size_t seq_len) {
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
void amplicon_missing_costs(double max_errors, size_t len, long k,
    long *skip, long *edge) {
  size_t max_t = len - min_truncated_len(len);
  for (size_t q = 0; q <= len; q++) {
    skip[q] = -1;
    edge[q] = -1;
  }
  for (size_t t = 0; t <= max_t; t++) {
    long k_t = k;
    if (t > 0) {
      k_t = amplicon_max_errors_of(max_errors, len - t) - 1;
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

/*
   Whether <hit> is better than <best> as the primer of an end of a read:
   the one which reaches furthest out (the one starting first at the
   left end, ending last at the right end), as the primer which started
   the amplicon is the outermost one and a primer-like sequence further
   inside is template; then the one with fewer errors, then the longer
   one. For a match inside the read (<end>: false), the one with fewer
   errors, then the longer one.

   The comparison is of the outer edge of the match, not of its inner
   edge: the primers of a scheme need not be equally long (those of
   influenza B are 20 to 28 bases), and the inner edge would prefer the
   shortest primer at an end to the one which really started the
   amplicon, however many errors that costs.
*/
static bool better_hit(const primer_hit_t *hit, const primer_hit_t *best,
    bool end, bool left) {
  if (end && left && hit->start != best->start)
    return hit->start < best->start;
  if (end && !left && hit->end != best->end)
    return hit->end > best->end;
  if (hit->n_errors != best->n_errors)
    return hit->n_errors < best->n_errors;
  return hit->end - hit->start > best->end - best->start;
}

/*
   Search the primer <p> at an end of <seq> (see primer_search_end()),
   with the costs <skip> and <edge> of its missing outer bases (see
   amplicon_missing_costs()) and the errors <k> allowed in the whole
   primer.

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
  out->keep_terminus = false;
  out->ref.mapped = false;
  bool end = partial_left || partial_right;
  size_t stride = self->max_primer_len + 1;
  for (size_t i = 0; i < n; i++) {
    if (!virus_set_has(self->types, primers[i].virus))
      continue;
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
void amplicon_search_layout(amplicon_analyzer_t *self, const char *seq,
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
  if (out->left.found && out->left.hit.start > left_window) {
    out->left.found = false;
    out->left.keep_terminus = false;
  }
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
  out->right.keep_terminus = false;
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

/* collect the matches of the primers of the scheme which lie between
   the primers of the two ends (see internal-matches.h) */
void amplicon_search_internal(amplicon_analyzer_t *self, const char *seq,
    size_t len, amplicon_t *out) {
  long lo = out->left.found ? out->left.hit.end + 1 : 0;
  long hi = out->right.found ? out->right.hit.start : (long)len;
  out->n_internal = internal_matches_find(self->internal, seq, lo, hi,
                                          out->internal);
}
