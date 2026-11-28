/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "edlib.h"
#include "edlib-iupac.h"
#include "alloc.h"
#include "numutils.h"
#include "primer-search.h"
#include "primers.h"

/* the character the read is padded with, to allow a primer to overhang
   it; it matches any base (it is the IUPAC code for "any base") */
#define PAD_CHAR 'N'

/* the cost of a start of the alignment which is not allowed, in the
   search at an end (see primer_search_end()) */
#define NO_START (INT_MAX / 2)

struct primer_search_t {
  char *padded;
  size_t padded_cap;
  /* the search at an end: whether two characters are equal (see
     IUPAC matching rules), the matrix of the dynamic programming, the
     primer and the read reversed (for the right end), the operations of
     the alignment path */
  bool eq[256][256];
  int *dp;
  size_t dp_cap;
  char *rev;
  size_t rev_cap;
  char *ops;
  size_t ops_cap;
};

primer_search_t* primer_search_new(void) {
  primer_search_t *self;
  XMALLOC(self, sizeof(primer_search_t));
  self->padded = NULL;
  self->padded_cap = 0;
  memset(self->eq, 0, sizeof(self->eq));
  for (int c = 0; c < 256; c++)
    self->eq[c][c] = true;
  for (size_t i = 0; i < n_iupac_equalities; i++) {
    unsigned char a = (unsigned char)iupac_equalities[i].first,
                  b = (unsigned char)iupac_equalities[i].second;
    self->eq[a][b] = true;
    self->eq[b][a] = true;
  }
  self->dp = NULL;
  self->dp_cap = 0;
  self->rev = NULL;
  self->rev_cap = 0;
  self->ops = NULL;
  self->ops_cap = 0;
  return self;
}

void primer_search_delete(primer_search_t *self) {
  if (self == NULL)
    return;
  free(self->padded);
  free(self->dp);
  free(self->rev);
  free(self->ops);
  free(self);
}

/* edit distance between the primer and target[start..end] */
static long distance(const char *primer, size_t primer_len,
    const char *target, long start, long end) {
  int n_pairs;
  const EdlibEqualityPair *pairs = edlib_iupac_equalities(&n_pairs);
  EdlibAlignResult result = edlibAlign(primer, (int)primer_len,
      target + start, (int)(end - start + 1),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE,
                          pairs, n_pairs));
  assert(result.status == EDLIB_STATUS_OK);
  long d = result.editDistance;
  edlibFreeAlignResult(result);
  return d;
}

/*
   Positions in <target> of the bases of the primer at positions
   <pos_a> <= <pos_b>, for the alignment of the whole primer to
   target[start..end] (the alignment path is computed only for a primer
   which was found).

   A primer base aligned to a base of the target gives its position; a
   primer base which is deleted in the target gives the position of the
   target base which follows the deletion.

   The operations of the alignment path of edlib consume, in the
   orientation used here (the primer is the query, the read the target):
   EDLIB_EDOP_MATCH and EDLIB_EDOP_MISMATCH a base of each,
   EDLIB_EDOP_INSERT a base of the primer only, EDLIB_EDOP_DELETE a base
   of the read only (the names are those of edlib, which calls an
   insertion what is seen as a deletion from the read; the unit tests
   below check the resulting positions on reads with an insertion and
   with a deletion).
*/
static void map_primer_positions(const char *primer, size_t primer_len,
    const char *target, long start, long end, size_t pos_a, size_t pos_b,
    long *out_a, long *out_b) {
  assert(pos_a <= pos_b);
  assert(pos_b <= primer_len);
  *out_a = start;
  *out_b = end + 1;
  if (pos_a == 0 && pos_b == 0)
    return;
  int n_pairs;
  const EdlibEqualityPair *pairs = edlib_iupac_equalities(&n_pairs);
  EdlibAlignResult result = edlibAlign(primer, (int)primer_len,
      target + start, (int)(end - start + 1),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH,
                          pairs, n_pairs));
  assert(result.status == EDLIB_STATUS_OK);
  long qpos = 0, tpos = start;
  bool got_a = pos_a == 0, got_b = false;
  for (int i = 0; i < result.alignmentLength && !got_b; i++) {
    /* a base of the read inserted at a sought position belongs to the
       part of the primer before it (it is consumed before the position
       is reached): an inserted base at the end of the tail is tail, and
       is thus removed, rather than kept as a base of the terminus */
    if (result.alignment[i] == EDLIB_EDOP_DELETE) {
      tpos++;
      continue;
    }
    if (!got_a && (size_t)qpos == pos_a) {
      *out_a = tpos;
      got_a = true;
    }
    if ((size_t)qpos == pos_b) {
      *out_b = tpos;
      got_b = true;
      break;
    }
    switch (result.alignment[i]) {
      case EDLIB_EDOP_MATCH:
      case EDLIB_EDOP_MISMATCH:
        qpos++;
        tpos++;
        break;
      case EDLIB_EDOP_INSERT:    /* a base of the primer, not of the read */
        qpos++;
        break;
      default:
        assert(false);
    }
  }
  /* a position at the end of the primer is reached after the last
     operation */
  if (!got_a && (size_t)qpos == pos_a)
    *out_a = tpos;
  if (!got_b && (size_t)qpos == pos_b)
    *out_b = tpos;
  edlibFreeAlignResult(result);
}

/*
   Search the primer in <target> (the read, of length <read_len>, padded
   by <pad_left> characters on the left and <target_len> - <read_len> -
   <pad_left> on the right); returns false if not found; otherwise the
   positions in the read (clipped to the read), whether the primer
   overhangs the read, the number of errors and the number of primer
   positions inside the read (in <in_read>, unless NULL)
*/
static bool search(const char *primer, size_t primer_len,
    size_t terminus_from, size_t terminus_len, const char *target,
    size_t target_len, long pad_left, long read_len, long max_errors,
    primer_hit_t *hit, long *in_read) {
  if (target_len == 0)
    return false;
  int n_pairs;
  const EdlibEqualityPair *pairs = edlib_iupac_equalities(&n_pairs);
  EdlibAlignResult result = edlibAlign(primer, (int)primer_len,
      target, (int)target_len,
      edlibNewAlignConfig((int)max_errors, EDLIB_MODE_HW, EDLIB_TASK_LOC,
                          pairs, n_pairs));
  assert(result.status == EDLIB_STATUS_OK);
  long n_errors = result.editDistance;
  if (n_errors < 0) {
    edlibFreeAlignResult(result);
    return false;
  }
  /* among equally good locations, prefer the (first) longest one */
  long start = result.startLocations[0], end = result.endLocations[0];
  for (int i = 1; i < result.numLocations; i++)
    if (result.endLocations[i] - result.startLocations[i] > end - start) {
      start = result.startLocations[i];
      end = result.endLocations[i];
    }
  edlibFreeAlignResult(result);
  /* extend the aligned region while the edit distance does not increase */
  while (start > 0 &&
         distance(primer, primer_len, target, start - 1, end) <= n_errors)
    start--;
  while (end < (long)target_len - 1 &&
         distance(primer, primer_len, target, start, end + 1) <= n_errors)
    end++;
  /* the bounds of the conserved terminus, before clipping the alignment
     to the read */
  long terminus_start, terminus_end;
  map_primer_positions(primer, primer_len, target, start, end,
                       terminus_from, terminus_from + terminus_len,
                       &terminus_start, &terminus_end);
  hit->partial = start < pad_left || end >= pad_left + read_len;
  start -= pad_left;
  end -= pad_left;
  long clipped_start = start > 0 ? start : 0;
  long clipped_end = end < read_len - 1 ? end : read_len - 1;
  if (in_read != NULL)
    *in_read = clipped_end - clipped_start + 1;
  hit->start = clipped_start;
  hit->end = clipped_end;
  hit->n_errors = n_errors;
  hit->terminus_start = terminus_start - pad_left;
  hit->terminus_end = terminus_end - pad_left;
  return true;
}

static void not_found(primer_hit_t *hit) {
  hit->start = -1;
  hit->end = -1;
  hit->partial = false;
  hit->n_errors = -1;
  hit->terminus_start = -1;
  hit->terminus_end = -1;
}

/* the search of a sequence whose last <terminus_len> bases are the
   conserved terminus (0: none) */
static bool run(primer_search_t *self, const char *seq, size_t seq_len,
    size_t terminus_from, size_t terminus_len, const char *read,
    size_t read_len, long max_errors, bool partial_left,
    bool partial_right, primer_hit_t *hit) {
  assert(self != NULL);
  assert(seq != NULL);
  assert(seq_len > 0);
  assert(terminus_from + terminus_len <= seq_len);
  assert(read != NULL);
  assert(hit != NULL);
  /* (unpadded: the whole alignment is in the read, in_read not needed) */
  if (search(seq, seq_len, terminus_from, terminus_len, read, read_len, 0,
             (long)read_len, max_errors, hit, NULL))
    return true;
  if (!partial_left && !partial_right) {
    not_found(hit);
    return false;
  }
  size_t max_overhang = seq_len / 2;
  size_t pad_left = partial_left ? max_overhang : 0;
  size_t pad_right = partial_right ? max_overhang : 0;
  size_t padded_len = pad_left + read_len + pad_right;
  XCHECKSPACEINC(self->padded, &self->padded_cap, padded_len + 1, 1024);
  memset(self->padded, PAD_CHAR, pad_left);
  memcpy(self->padded + pad_left, read, read_len);
  memset(self->padded + pad_left + read_len, PAD_CHAR, pad_right);
  self->padded[padded_len] = '\0';
  long in_read;
  if (!search(seq, seq_len, terminus_from, terminus_len, self->padded,
              padded_len, (long)pad_left, (long)read_len, max_errors, hit,
              &in_read) ||
      hit->n_errors > floordiv(max_errors * in_read, (long)seq_len)) {
    not_found(hit);
    return false;
  }
  return true;
}

bool primer_search_run(primer_search_t *self, const primer_t *primer,
    const char *read, size_t read_len, long max_errors, bool partial_left,
    bool partial_right, primer_hit_t *hit) {
  assert(primer != NULL);
  return run(self, primer->seq, primer->seq_len, primer->terminus_from,
             primer->terminus_len, read, read_len, max_errors,
             partial_left, partial_right, hit);
}

bool primer_search_run_seq(primer_search_t *self, const char *seq,
    size_t seq_len, const char *read, size_t read_len, long max_errors,
    bool partial_left, bool partial_right, primer_hit_t *hit) {
  /* no terminus: the whole alignment is tail (see primer-search.h) */
  return run(self, seq, seq_len, seq_len, 0, read, read_len, max_errors,
             partial_left, partial_right, hit);
}

/* the cost of a start of the alignment after <q> missing bases, from a
   cost array of primer_search_end() */
static int start_cost(const long *cost, size_t q, long max_cost) {
  return cost[q] < 0 || cost[q] > max_cost ? NO_START : (int)cost[q];
}

/*
   The position in the read of the primer base <pos>, from the
   operations <ops> (<n_ops>, in the order of the read, or in the
   opposite order if <backward>) of an alignment which begins with the
   primer base <qpos> at the read position <tpos> (<pos> is at least
   <qpos> and at most the position after the last aligned base); as in
   map_primer_positions(), a read base inserted before a primer base is
   counted before it
*/
static long map_position(const char *ops, size_t n_ops, bool backward,
    long qpos, long tpos, long pos) {
  for (size_t i = 0; i < n_ops; i++) {
    char op = backward ? ops[n_ops - 1 - i] : ops[i];
    if (op == EDLIB_EDOP_DELETE) {
      tpos++;
      continue;
    }
    if (qpos == pos)
      return tpos;
    qpos++;
    if (op != EDLIB_EDOP_INSERT)
      tpos++;
  }
  return tpos;
}

/*
   The path of the alignment ending with the read base <end> - 1 (see
   primer_search_end(): the primer <a> (length <m>) and the read <b>
   from the outer side, the matrix in self->dp, of width <w>), from its
   inner end, into self->ops; the preferred operation is a match, then
   the start of the alignment (fewer errors), then a mismatch, then a
   gap. Returns the number of operations; <out_q> and <out_j> are set to
   the first primer and read positions of the alignment, <out_cost> to
   the cost of its start.
*/
static size_t trace_back(primer_search_t *self, const char *a,
    const char *b, size_t m, size_t w, const long *skip_cost,
    const long *edge_cost, long max_cost, size_t end, size_t *out_q,
    size_t *out_j, int *out_cost) {
  const int *d = self->dp;
  size_t n_ops = 0, q = m, j = end;
  int cost;
  for (;;) {
    int v = d[q * w + j];
    const int *diag = q > 0 && j > 0 ? d + (q - 1) * w + j - 1 : NULL;
    bool eq = diag != NULL &&
              self->eq[(unsigned char)a[q - 1]][(unsigned char)b[j - 1]];
    if (eq && *diag == v) {
      self->ops[n_ops++] = EDLIB_EDOP_MATCH;
      q--;
      j--;
      continue;
    }
    cost = q == 0 ? 0 : start_cost(j == 0 ? edge_cost : skip_cost, q,
                                   max_cost);
    if (cost == v)
      break;
    if (diag != NULL && !eq && *diag + 1 == v) {
      self->ops[n_ops++] = EDLIB_EDOP_MISMATCH;
      q--;
      j--;
    } else if (q > 0 && d[(q - 1) * w + j] + 1 == v) {
      // a primer base deleted in the read
      self->ops[n_ops++] = EDLIB_EDOP_INSERT;
      q--;
    } else {
      // a read base inserted
      assert(j > 0 && d[q * w + j - 1] + 1 == v);
      self->ops[n_ops++] = EDLIB_EDOP_DELETE;
      j--;
    }
  }
  *out_q = q;
  *out_j = j;
  *out_cost = cost;
  return n_ops;
}

bool primer_search_end(primer_search_t *self, const primer_t *primer,
    const char *read, size_t read_len, bool left, const long *skip_cost,
    const long *edge_cost, long max_cost, primer_hit_t *hit) {
  assert(self != NULL);
  assert(primer != NULL);
  assert(primer->seq_len > 0);
  assert(read != NULL);
  assert(skip_cost != NULL && skip_cost[0] == 0);
  assert(edge_cost != NULL && edge_cost[0] == 0);
  assert(hit != NULL);
  size_t m = primer->seq_len, n = read_len;
  if (n == 0 || max_cost < 0) {
    not_found(hit);
    return false;
  }
  // the primer (a) and the read (b) from the outer side
  const char *a = primer->seq, *b = read;
  if (!left) {
    XCHECKSPACEINC(self->rev, &self->rev_cap, m + n, 256);
    for (size_t i = 0; i < m; i++)
      self->rev[i] = primer->seq[m - 1 - i];
    for (size_t i = 0; i < n; i++)
      self->rev[m + i] = read[n - 1 - i];
    a = self->rev;
    b = self->rev + m;
  }
  /*
     d[q * w + j]: the lowest cost of an alignment of the first q primer
     bases (from the outer side, some of them missing at a cost) ending
     with the read base j - 1; the read bases before the alignment are
     free
  */
  size_t w = n + 1;
  XCHECKSPACEINC(self->dp, &self->dp_cap, sizeof(int) * (m + 1) * w,
                 sizeof(int) * 4096);
  int *d = self->dp;
  for (size_t j = 0; j <= n; j++)
    d[j] = 0;
  for (size_t q = 1; q <= m; q++) {
    int *row = d + q * w, *up = row - w;
    const bool *eq = self->eq[(unsigned char)a[q - 1]];
    int v = start_cost(edge_cost, q, max_cost);
    row[0] = up[0] + 1 < v ? up[0] + 1 : v;
    int skip = start_cost(skip_cost, q, max_cost);
    for (size_t j = 1; j <= n; j++) {
      v = up[j - 1] + !eq[(unsigned char)b[j - 1]];
      if (up[j] + 1 < v)
        v = up[j] + 1;
      if (row[j - 1] + 1 < v)
        v = row[j - 1] + 1;
      if (skip < v)
        v = skip;
      row[j] = v;
    }
  }
  /*
     the end with the lowest cost; among equal ones, the outermost. Two
     ends of equal cost differ by read bases which the read alone cannot
     assign: the inner one takes them as primer, the outer one as
     template. The outer one is taken, so that in doubt a base of
     template is kept rather than removed: e.g. where the last primer
     base, G, meets a template G and the read has an error in the primer,
     an alignment ending one base further in costs the same, and would
     remove the template G (the price: a primer base with an error in
     the read, e.g. a substituted last base, is kept instead)
  */
  const int *last = d + m * w;
  size_t best = 0;
  for (size_t j = 1; j <= n; j++)
    if (last[j] <= max_cost && (best == 0 || last[j] < last[best]))
      best = j;
  if (best == 0) {
    not_found(hit);
    return false;
  }
  XCHECKSPACEINC(self->ops, &self->ops_cap, m + n, 256);
  size_t q, j;
  int cost;
  size_t n_ops = trace_back(self, a, b, m, w, skip_cost, edge_cost,
                            max_cost, best, &q, &j, &cost);
  /* the alignment: primer bases q to m - 1 (from the outer side) with
     the read bases j to best - 1 (from the outer side) */
  bool edge = j == 0 && q > 0;
  hit->partial = edge;
  hit->n_errors = last[best] - cost;
  long skipped = (long)q, len = (long)n;
  long pa = (long)primer->terminus_from,
       pb = (long)(primer->terminus_from + primer->terminus_len);
  if (left) {
    hit->start = (long)j;
    hit->end = (long)best - 1;
    // the operations are in the opposite order of the read
    hit->terminus_start = pa < skipped ?
        (edge ? pa - skipped : hit->start) :
        map_position(self->ops, n_ops, true, skipped, hit->start, pa);
    hit->terminus_end = pb < skipped ?
        (edge ? pb - skipped : hit->start) :
        map_position(self->ops, n_ops, true, skipped, hit->start, pb);
  } else {
    hit->start = len - (long)best;
    hit->end = len - 1 - (long)j;
    // the last aligned primer base is m - 1 - skipped
    long after = (long)m - skipped;
    hit->terminus_start = pa > after ?
        (edge ? hit->end + 1 + pa - after : hit->end + 1) :
        map_position(self->ops, n_ops, false, 0, hit->start, pa);
    hit->terminus_end = pb > after ?
        (edge ? hit->end + 1 + pb - after : hit->end + 1) :
        map_position(self->ops, n_ops, false, 0, hit->start, pb);
  }
  return true;
}
