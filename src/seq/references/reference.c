/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "fasta.h"
#include "kmer-index.h"
#include "pattern-scan.h"
#include "primers.h"
#include "reference.h"
#include "reference-index.h"
#include "segments.h"

struct reference_set_t {
  const reference_index_t *index;
  const fasta_t *fasta;
  size_t probe_len;
  double max_error_rate;
  bool exhaustive;
  bool seeds_only;
  /* if not NULL, whether each reference is considered (see
     reference_set_restrict) */
  bool *allowed;
  /* if > 0, only the references of this segment are considered (see
     reference_set_restrict_segment); 0: all */
  int restrict_segment;
  /* the reverse complement of the probe, and the probe reversed
     (probe_len + 1 bytes each) */
  char *rc, *reversed;
  /* the probe prepared for its search in a window of a reference
     (pattern-scan.h): on the strand of the reference, read backwards
     (the probe reversed), and on the other strand (its reverse
     complement) */
  pattern_scan_t *scan_forward, *scan_reverse;
  /* the candidate diagonals (a buffer of hits_size bytes) */
  kmer_index_hit_t *hits;
  size_t hits_size;
};

reference_set_t* reference_set_new(const reference_index_t *index,
    size_t probe_len, double max_error_rate) {
  assert(index != NULL);
  assert(probe_len >= REFERENCE_MIN_PROBE_LEN);
  assert(max_error_rate >= 0 && max_error_rate < 1);
  reference_set_t *self;
  XMALLOC(self, sizeof(reference_set_t));
  self->index = index;
  self->fasta = reference_index_fasta(index);
  self->probe_len = probe_len;
  self->max_error_rate = max_error_rate;
  self->exhaustive = false;
  self->seeds_only = false;
  self->allowed = NULL;
  self->restrict_segment = 0;
  XMALLOC(self->rc, probe_len + 1);
  XMALLOC(self->reversed, probe_len + 1);
  self->scan_forward = pattern_scan_new("N", 1);
  self->scan_reverse = pattern_scan_new("N", 1);
  self->hits = NULL;
  self->hits_size = 0;
  return self;
}

void reference_set_set_seeds_only(reference_set_t *self, bool seeds_only) {
  self->seeds_only = seeds_only;
}

void reference_set_set_exhaustive(reference_set_t *self, bool exhaustive) {
  self->exhaustive = exhaustive;
}

void reference_set_restrict(reference_set_t *self, const bool *allowed) {
  free(self->allowed);
  self->allowed = NULL;
  size_t n = self->fasta->n_records;
  if (allowed == NULL || n == 0)
    return;
  XMALLOC(self->allowed, sizeof(bool) * n);
  memcpy(self->allowed, allowed, sizeof(bool) * n);
}

void reference_set_restrict_segment(reference_set_t *self, int segment) {
  assert(segment >= 0 && segment <= SEGMENTS_N);
  self->restrict_segment = segment;
}

/* whether the reference <ref> is considered, given the restrictions to
   some references and to a segment */
static bool ref_allowed(const reference_set_t *self, size_t ref) {
  return (self->allowed == NULL || self->allowed[ref]) &&
         (self->restrict_segment == 0 ||
          reference_index_segment(self->index, ref) ==
          self->restrict_segment);
}

bool reference_set_allows(const reference_set_t *self, size_t ref) {
  assert(self != NULL);
  return ref_allowed(self, ref);
}

void reference_set_delete(reference_set_t *self) {
  if (self == NULL)
    return;
  free(self->rc);
  free(self->reversed);
  pattern_scan_delete(self->scan_forward);
  pattern_scan_delete(self->scan_reverse);
  free(self->hits);
  free(self->allowed);
  free(self);
}

size_t reference_set_size(const reference_set_t *self) {
  return self->fasta->n_records;
}

size_t reference_set_probe_len(const reference_set_t *self) {
  return self->probe_len;
}

const fasta_t* reference_set_fasta(const reference_set_t *self) {
  return self->fasta;
}

/* distance of a hit from the terminus from which the read would start */
static long terminus_distance(const reference_hit_t *hit) {
  return hit->reverse ? hit->ref_len - 1 - hit->ref_start : hit->ref_start;
}

/* prepare the search of <probe> (length <len>) on the strand of the
   references (see align_window()) */
static void set_forward(reference_set_t *self, const char *probe,
    size_t len) {
  for (size_t i = 0; i < len; i++)
    self->reversed[i] = probe[len - 1 - i];
  self->reversed[len] = '\0';
  pattern_scan_set(self->scan_forward, self->reversed, len);
}

/* the same on the other strand; the reverse complement of the probe is
   left in self->rc */
static void set_reverse(reference_set_t *self, const char *probe,
    size_t len) {
  primers_reverse_complement(probe, len, self->rc);
  self->rc[len] = '\0';
  pattern_scan_set(self->scan_reverse, self->rc, len);
}

/* whether two references are of the same segment and of the same
   subtype (false if either has no subtype in its description) */
static bool same_subtype(const reference_set_t *self, size_t a, size_t b) {
  int subtype = reference_index_subtype(self->index, a);
  return subtype >= 0 && subtype == reference_index_subtype(self->index, b)
         && reference_index_segment(self->index, a) ==
            reference_index_segment(self->index, b);
}

/* the match of the probe at <ref_start> of the reference <ref_index>
   (on the other strand if <reverse>) with <n_errors> errors: stored in
   <best> if better than it (see align_window()) */
static void consider_match(const reference_set_t *self, size_t ref_index,
    bool reverse, long ref_start, long n_errors, reference_hit_t *best) {
  reference_hit_t hit = {.mapped = true, .ref_index = ref_index,
      .ref_len = (long)self->fasta->records[ref_index].seq_len,
      .reverse = reverse, .ref_start = ref_start, .n_errors = n_errors,
      .ambiguous = false};
  if (!best->mapped || hit.n_errors < best->n_errors) {
    *best = hit;
  } else if (hit.n_errors == best->n_errors) {
    /* as good as the best one: nearer to the terminus, it replaces it;
       on a reference of another segment or subtype, the probe does not
       tell the two apart. Two references of one segment and subtype
       (the lineages of a subtype) are not told apart either, but the
       probe does say which subtype it belongs to, which is what the
       choice of the reference needs: it is not ambiguous, and the
       match kept stands for both. */
    bool ambiguous = best->ambiguous ||
        (hit.ref_index != best->ref_index &&
         !same_subtype(self, hit.ref_index, best->ref_index));
    if (terminus_distance(&hit) < terminus_distance(best))
      *best = hit;
    best->ambiguous = ambiguous;
  }
}

/*
   Align the probe prepared by set_forward() (or, for <reverse>, by
   set_reverse()) to the positions <from> to <to> of the reference
   <ref_index> with at most <max_errors> errors; if found, and better
   than <best> (fewer errors, or as many and nearer to the terminus, see
   reference.h), store the match in <best>. The position of a match is
   that of the first base of the probe: where the alignment begins or,
   for <reverse>, where it ends.

   (A run of N in a reference matches any probe: the references are
   expected to be complete sequences of the segments.)
*/
static void align_window(reference_set_t *self, size_t ref_index,
    long from, long to, bool reverse, long max_errors,
    reference_hit_t *best) {
  const fasta_record_t *ref = self->fasta->records + ref_index;
  if (from < 0)
    from = 0;
  if (to > (long)ref->seq_len - 1)
    to = (long)ref->seq_len - 1;
  if (to < from)
    return;
  /* the location nearest to the terminus (there can be several, e.g.
     for a repeat, or for equally good alignments of different lengths):
     the first start on the strand of the reference, the last end on the
     other one */
  long n_errors, first, last;
  if (!pattern_scan_best(reverse ? self->scan_reverse : self->scan_forward,
                         ref->seq + from, (size_t)(to - from + 1),
                         max_errors, !reverse, &n_errors, &first, &last))
    return;
  consider_match(self, ref_index, reverse, from + (reverse ? last : first),
                 n_errors, best);
}

/* align_window() of the probe prepared by set_forward() and by
   set_reverse() to the whole reference <ref_index>: first on its
   strand, then on the other one, in one pass over the reference */
static void align_reference(reference_set_t *self, size_t ref_index,
    long max_errors, reference_hit_t *best) {
  const fasta_record_t *ref = self->fasta->records + ref_index;
  pattern_scan_best_t forward, reverse;
  pattern_scan_best_both(self->scan_forward, self->scan_reverse, ref->seq,
                         ref->seq_len, max_errors, &forward, &reverse);
  /* (a match on the other strand with more errors than the one on the
     strand of the reference does not replace it) */
  if (forward.found)
    consider_match(self, ref_index, false, forward.first, forward.n_errors,
                   best);
  if (reverse.found)
    consider_match(self, ref_index, true, reverse.last, reverse.n_errors,
                   best);
}

/*
   The seeded search of <query> (one strand of the probe): its k-mers give
   candidate diagonals in the references (kmer-index.h, sorted by
   reference and diagonal); the diagonals of a reference closer than the
   probe length are merged into one window, and each window is aligned
   (with <max_errors> errors of slack on both sides)
*/
static void seeded_search(reference_set_t *self, const char *query,
    size_t len, bool reverse, long max_errors, reference_hit_t *best) {
  size_t n = kmer_index_diagonals(reference_index_kmers(self->index),
                                  query, len, REFERENCE_MAX_KMER_OCC,
                                  &self->hits, &self->hits_size);
  for (size_t i = 0; i < n;) {
    uint32_t ref = self->hits[i].seq;
    if (!ref_allowed(self, ref)) {
      size_t j = i + 1;
      while (j < n && self->hits[j].seq == ref)
        j++;
      i = j;
      continue;
    }
    long first = self->hits[i].diagonal;
    long last = first;
    size_t j = i + 1;
    while (j < n && self->hits[j].seq == ref &&
           self->hits[j].diagonal <= last + (long)len) {
      last = self->hits[j].diagonal;
      j++;
    }
    // a worse match than the best one so far is not needed (an equally
    // good one is, for the choice of the location)
    long k = best->mapped ? best->n_errors : max_errors;
    align_window(self, ref, first - max_errors,
                 last + (long)len - 1 + max_errors, reverse, k, best);
    i = j;
  }
}

/*
   Whether the first <len> bases of <probe> can be located: at least
   REFERENCE_MIN_PROBE_LEN of them other than N (an N matches every base
   of a reference, so that a probe of Ns would be located, with no
   error, wherever the search looks first)
*/
static bool locatable(const char *probe, size_t len) {
  size_t known = 0;
  for (size_t i = 0; i < len; i++)
    if (probe[i] != 'N' && probe[i] != 'n')
      known++;
  return known >= REFERENCE_MIN_PROBE_LEN;
}

/* locate a probe (see reference_set_locate()); <exhaustive>: without the
   index; <seeds_only>: without the exhaustive search of a probe which the
   index does not locate */
static void locate(reference_set_t *self, const char *probe,
    size_t probe_len, bool exhaustive, bool seeds_only,
    reference_hit_t *hit) {
  assert(probe != NULL);
  assert(hit != NULL);
  memset(hit, 0, sizeof(reference_hit_t));
  hit->mapped = false;
  hit->ref_start = -1;
  hit->n_errors = -1;
  size_t len = probe_len < self->probe_len ? probe_len : self->probe_len;
  if (!locatable(probe, len))
    return;
  set_forward(self, probe, len);
  set_reverse(self, probe, len);
  long max_errors = (long)((double)len * self->max_error_rate);
  if (!exhaustive) {
    seeded_search(self, probe, len, false, max_errors, hit);
    seeded_search(self, self->rc, len, true, max_errors, hit);
    // nothing found from the k-mers (e.g. an error in each of them):
    // the exhaustive search follows, unless seeds only
    if (hit->mapped || seeds_only)
      return;
  }
  for (size_t i = 0; i < self->fasta->n_records; i++) {
    if (!ref_allowed(self, i))
      continue;
    // a worse match than the best one so far is not needed
    align_reference(self, i, hit->mapped ? hit->n_errors : max_errors, hit);
  }
}

void reference_set_locate(reference_set_t *self, const char *probe,
    size_t probe_len, reference_hit_t *hit) {
  locate(self, probe, probe_len, self->exhaustive, self->seeds_only, hit);
}

void reference_set_locate_seeded(reference_set_t *self, const char *probe,
    size_t probe_len, reference_hit_t *hit) {
  locate(self, probe, probe_len, false, true, hit);
}

void reference_set_locate_within(reference_set_t *self, const char *probe,
    size_t probe_len, size_t ref_index, bool reverse, long from, long to,
    reference_hit_t *hit) {
  assert(probe != NULL);
  assert(hit != NULL);
  assert(ref_index < self->fasta->n_records);
  memset(hit, 0, sizeof(reference_hit_t));
  hit->mapped = false;
  hit->ref_start = -1;
  hit->n_errors = -1;
  size_t len = probe_len < self->probe_len ? probe_len : self->probe_len;
  if (!locatable(probe, len))
    return;
  if (reverse)
    set_reverse(self, probe, len);
  else
    set_forward(self, probe, len);
  align_window(self, ref_index, from, to, reverse,
               (long)((double)len * self->max_error_rate), hit);
}

void reference_set_locate_in(reference_set_t *self, const char *probe,
    size_t probe_len, size_t ref_index, bool reverse,
    reference_hit_t *hit) {
  assert(ref_index < self->fasta->n_records);
  reference_set_locate_within(self, probe, probe_len, ref_index, reverse,
      0, (long)self->fasta->records[ref_index].seq_len - 1, hit);
}

bool reference_hit_is_terminal(const reference_hit_t *hit,
    long terminal_window) {
  assert(hit->mapped);
  return terminus_distance(hit) < terminal_window;
}
