/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "anchor-site.h"
#include "codons.h"
#include "examples-table.h"
#include "alloc.h"
#include "edlib.h"
#include "edlib-iupac.h"
#include "iupac.h"
#include "sequtil.h"
#include "test.h"

static char upper(char c) {
  return c >= 'a' && c <= 'z' ? (char)(c - 'a' + 'A') : c;
}

/* whether the <avail> residues at <aa> begin with <peptide>, in which
   a '.' matches any residue (as codons_find_peptide()) */
static bool peptide_matches(const char *peptide, const char *aa,
    size_t avail) {
  size_t i = 0;
  for (const char *p = peptide; *p != '\0'; p++, i++) {
    if (i >= avail)
      return false;
    if (*p != '.' && *p != aa[i])
      return false;
  }
  return true;
}

size_t anchor_peptide_occurrences(const char *seq, size_t seq_len,
    const char *peptide, long *first) {
  assert(peptide != NULL);
  assert(first != NULL);
  *first = -1;
  if (seq == NULL || seq_len < 3 || peptide[0] == '\0')
    return 0;
  char *protein;
  XMALLOC(protein, seq_len / 3 + 1);
  size_t n_found = 0;
  size_t plen = strlen(peptide);
  for (size_t frame = 0; frame < 3 && seq_len > frame; frame++) {
    size_t n = codons_translate(seq + frame, seq_len - frame, protein);
    for (size_t i = 0; i + plen <= n; i++)
      if (peptide_matches(peptide, protein + i, n - i)) {
        if (*first < 0)
          *first = (long)frame + 3 * (long)i;
        n_found++;
      }
  }
  free(protein);
  return n_found;
}

bool anchor_bases_alone(const anchor_t *self, const char *seq,
    size_t seq_len, long pos, double max_error_rate) {
  assert(self != NULL);
  if (!self->found || seq == NULL || pos < 0)
    return false;
  anchor_hit_t hit;
  size_t after = (size_t)pos + self->len;
  return after <= seq_len &&
         !anchor_find(self, seq, (size_t)pos, max_error_rate, &hit) &&
         !anchor_find(self, seq + after, seq_len - after, max_error_rate,
                      &hit);
}

bool anchor_of_peptide(anchor_t *self, size_t ref, const char *seq,
    size_t seq_len, const char *peptide, size_t len) {
  return anchor_at(self, ref, seq, seq_len,
                   codons_find_peptide(seq, seq_len, peptide), len);
}

bool anchor_at(anchor_t *self, size_t ref, const char *seq, size_t seq_len,
    long pos, size_t len) {
  assert(self != NULL);
  memset(self, 0, sizeof(anchor_t));
  if (len == 0 || len > ANCHOR_MAX_LEN)
    return false;
  if (pos < 0 || (size_t)pos + len > seq_len)
    return false;
  self->found = true;
  self->ref = ref;
  self->ref_pos = pos;
  self->len = len;
  for (size_t i = 0; i < len; i++)
    self->seq[i] = upper(seq[pos + i]);
  self->seq[len] = '\0';
  reverse_complement(self->seq, len, self->rev);
  return true;
}

/* the best infix alignment of the <query_len> characters of <query> in
   the <target_len> characters of <target>, with at most <max_errors>
   errors: its edit distance, or -1 if there is none; <start> and <end>
   are the 0-based first and last positions of the alignment in the
   target */
static long search(const char *query, size_t query_len, const char *target,
    size_t target_len, long max_errors, long *start, long *end) {
  if (target_len == 0 || target_len < query_len)
    return -1;
  int n_pairs;
  const EdlibEqualityPair *pairs = edlib_iupac_equalities(&n_pairs);
  EdlibAlignResult r = edlibAlign(query, (int)query_len, target,
      (int)target_len,
      edlibNewAlignConfig((int)max_errors, EDLIB_MODE_HW, EDLIB_TASK_LOC,
                          pairs, n_pairs));
  assert(r.status == EDLIB_STATUS_OK);
  long distance = r.editDistance;
  if (distance < 0) {
    edlibFreeAlignResult(r);
    return -1;
  }
  /* among equally good locations, the (first) longest one, as in the
     search of the primers and of the adaptors */
  *start = r.startLocations[0];
  *end = r.endLocations[0];
  for (int i = 1; i < r.numLocations; i++)
    if (r.endLocations[i] - r.startLocations[i] > *end - *start) {
      *start = r.startLocations[i];
      *end = r.endLocations[i];
    }
  edlibFreeAlignResult(r);
  return distance;
}

/* the stretch of a read of <read_len> bases from <from> to <to> (0-based,
   the last one included), clamped to the read: its first position in
   <first>, and its length (0 if nothing of it is in the read) */
static size_t clamp(long from, long to, size_t read_len, long *first) {
  if (from < 0)
    from = 0;
  if (to > (long)read_len - 1)
    to = (long)read_len - 1;
  *first = from;
  return to >= from ? (size_t)(to - from + 1) : 0;
}

bool anchor_find_within(const anchor_t *self, const char *read,
    size_t read_len, long fwd_from, long fwd_to, long rev_from, long rev_to,
    double max_error_rate, anchor_hit_t *hit) {
  assert(self != NULL);
  assert(hit != NULL);
  memset(hit, 0, sizeof(anchor_hit_t));
  if (!self->found || read == NULL || read_len == 0)
    return false;
  long fwd_first, rev_first;
  size_t fwd_len = clamp(fwd_from, fwd_to, read_len, &fwd_first);
  size_t rev_len = clamp(rev_from, rev_to, read_len, &rev_first);
  long max_errors = (long)(max_error_rate * (double)self->len);
  long fwd_start = 0, fwd_end = 0, rev_start = 0, rev_end = 0;
  long fwd = search(self->seq, self->len, read + fwd_first, fwd_len,
                    max_errors, &fwd_start, &fwd_end);
  long rev = search(self->rev, self->len, read + rev_first, rev_len,
                    max_errors, &rev_start, &rev_end);
  if (fwd < 0 && rev < 0)
    return false;
  /* as good in both orientations: which side of the match a window is
     on is not decided */
  if (fwd >= 0 && rev >= 0 && fwd == rev)
    return false;
  hit->found = true;
  hit->reverse = !(rev < 0 || (fwd >= 0 && fwd < rev));
  hit->start = hit->reverse ? rev_first + rev_start : fwd_first + fwd_start;
  hit->end = hit->reverse ? rev_first + rev_end : fwd_first + fwd_end;
  return true;
}

bool anchor_find(const anchor_t *self, const char *read, size_t read_len,
    double max_error_rate, anchor_hit_t *hit) {
  return anchor_find_within(self, read, read_len, 0, (long)read_len - 1, 0,
                            (long)read_len - 1, max_error_rate, hit);
}

bool anchor_window_span(const anchor_hit_t *hit, size_t read_len,
    long offset, size_t n, long *from, long *to) {
  assert(hit != NULL);
  assert(from != NULL);
  assert(to != NULL);
  if (!hit->found || n == 0)
    return false;
  /* the read runs as the reference: the window begins <offset> bases
     after the first base of the match. The read the other way round:
     the first base of the anchor there is the last base of the match
     here, and the window runs backwards from it */
  long at = hit->reverse ? hit->end - offset - (long)n + 1
                         : hit->start + offset;
  if (at < 0 || at + (long)n > (long)read_len)
    return false;
  *from = at;
  *to = at + (long)n - 1;
  return true;
}

bool anchor_bases_window(const anchor_hit_t *hit, const char *read,
    size_t read_len, long offset, size_t n, char *out) {
  assert(out != NULL);
  long from, to;
  if (read == NULL ||
      !anchor_window_span(hit, read_len, offset, n, &from, &to))
    return false;
  if (hit->reverse) {
    reverse_complement(read + from, n, out);
    return true;
  }
  for (size_t i = 0; i < n; i++)
    out[i] = upper(read[from + (long)i]);
  out[n] = '\0';
  return true;
}

bool anchor_bases_before(const anchor_hit_t *hit, const char *read,
    size_t read_len, size_t n, char *out) {
  return anchor_bases_window(hit, read, read_len, -(long)n, n, out);
}

bool anchor_bases_at(const anchor_hit_t *hit, const char *read,
    size_t read_len, size_t n, char *out) {
  return anchor_bases_window(hit, read, read_len, 0, n, out);
}

bool anchor_bases_after(const anchor_hit_t *hit, const char *read,
    size_t read_len, size_t skip, size_t n, char *out) {
  if (!hit->found)
    return false;
  return anchor_bases_window(hit, read, read_len,
      (long)anchor_hit_len(hit) + (long)skip, n, out);
}

/* --- unit tests --- */

/* the fusion peptide which begins HA2, and the cleavage site before it
   in a highly pathogenic H5 (metadata/examples.yaml) */
#define FUSION_NT EXAMPLE_FUSION_PEPTIDE
#define SITE_NT EXAMPLE_CLEAVAGE_SITE_H5_HPAI

void anchor_site_test(void) {
  const char *ha = "ATGGAAAAA" SITE_NT FUSION_NT "TTTATAGAGGGA";
  anchor_t a;
  test_assert(anchor_of_peptide(&a, 3, ha, strlen(ha), "G..GAIAG", 24));
  test_assert(a.found);
  test_eq_size(a.ref, 3);
  test_eq_long(a.ref_pos, 9 + (long)strlen(SITE_NT));
  test_eq_str(a.seq, FUSION_NT);
  test_eq_str(a.rev, "GCCTGCTATAGCCCCAAACAGTCC");
  /* a peptide which is not there, and an anchor past the end */
  anchor_t b;
  test_assert(!anchor_of_peptide(&b, 0, ha, strlen(ha), "WWWWW", 24));
  test_assert(!anchor_of_peptide(&b, 0, ha, strlen(ha), "G..GAIAG", 48));
  test_assert(!anchor_of_peptide(&b, 0, ha, strlen(ha), "G..GAIAG", 0));

  char out[32];
  anchor_hit_t hit;
  /* the read as the reference: the last 9 codons before HA2 */
  test_assert(anchor_find(&a, ha, strlen(ha), 0.2, &hit));
  test_assert(!hit.reverse);
  test_assert(anchor_bases_before(&hit, ha, strlen(ha), 27, out));
  test_eq_str(out, SITE_NT);
  /* and the bases after it, the control windows: the first 6 of HA2
     after the anchor, then the next ones */
  test_assert(anchor_bases_after(&hit, ha, strlen(ha), 0, 6, out));
  test_eq_str(out, "TTTATA");
  test_assert(anchor_bases_after(&hit, ha, strlen(ha), 6, 6, out));
  test_eq_str(out, "GAGGGA");
  /* the same read, reverse complemented: the same windows come back */
  char rc[128];
  reverse_complement(ha, strlen(ha), rc);
  test_assert(anchor_find(&a, rc, strlen(rc), 0.2, &hit));
  test_assert(hit.reverse);
  test_assert(anchor_bases_before(&hit, rc, strlen(rc), 27, out));
  test_eq_str(out, SITE_NT);
  test_assert(anchor_bases_after(&hit, rc, strlen(rc), 0, 6, out));
  test_eq_str(out, "TTTATA");
  test_assert(anchor_bases_after(&hit, rc, strlen(rc), 6, 6, out));
  test_eq_str(out, "GAGGGA");
  /* searched in a part of the read only: the positions are those of the
     whole read; outside the part the anchor is not seen, and each
     orientation is looked for in its own part */
  long at = 9 + (long)strlen(SITE_NT);
  test_assert(anchor_find_within(&a, ha, strlen(ha), at - 5, at + 30, 0, 0,
                                 0.2, &hit));
  test_assert(!hit.reverse);
  test_eq_long(hit.start, at);
  test_eq_long(hit.end, at + 23);
  test_assert(!anchor_find_within(&a, ha, strlen(ha), 0, at - 1, 0,
                                  (long)strlen(ha) - 1, 0.2, &hit));
  test_assert(anchor_find_within(&a, rc, strlen(rc), 0, 0, -10, 1000, 0.2,
                                 &hit));
  test_assert(hit.reverse);
  test_assert(!anchor_find_within(&a, rc, strlen(rc), 0,
                                  (long)strlen(rc) - 1, 0, 10, 0.2, &hit));
  /* a read with errors in the anchor still finds it */
  char noisy[128];
  snprintf(noisy, sizeof(noisy), "ATGGAAAAA" SITE_NT
           "GGACTCTTTGGGGCTTTAGCAGGC" "TTTATAGAGGGA");
  test_assert(anchor_find(&a, noisy, strlen(noisy), 0.2, &hit));
  test_assert(anchor_bases_before(&hit, noisy, strlen(noisy), 27, out));
  test_eq_str(out, SITE_NT);
  /* a read which does not reach far enough on either side */
  const char *short_read = "AAAGA" FUSION_NT "TTT";
  test_assert(anchor_find(&a, short_read, strlen(short_read), 0.2, &hit));
  test_assert(!anchor_bases_before(&hit, short_read, strlen(short_read),
                                   27, out));
  test_assert(anchor_bases_before(&hit, short_read, strlen(short_read),
                                  5, out));
  test_eq_str(out, "AAAGA");
  test_assert(anchor_bases_after(&hit, short_read, strlen(short_read), 0,
                                 3, out));
  test_eq_str(out, "TTT");
  test_assert(!anchor_bases_after(&hit, short_read, strlen(short_read), 0,
                                  4, out));
  test_assert(!anchor_bases_after(&hit, short_read, strlen(short_read), 3,
                                  3, out));
  /* the anchor itself, as the read carries it */
  test_assert(anchor_find(&a, ha, strlen(ha), 0.2, &hit));
  test_assert(anchor_bases_at(&hit, ha, strlen(ha), 24, out));
  test_eq_str(out, FUSION_NT);
  test_assert(anchor_find(&a, rc, strlen(rc), 0.2, &hit));
  test_assert(anchor_bases_at(&hit, rc, strlen(rc), 24, out));
  test_eq_str(out, FUSION_NT);
  /* a read whose anchor carries an error shows that error */
  test_assert(anchor_find(&a, noisy, strlen(noisy), 0.2, &hit));
  test_assert(anchor_bases_at(&hit, noisy, strlen(noisy), 24, out));
  test_eq_str(out, "GGACTCTTTGGGGCTTTAGCAGGC");
  test_assert(anchor_find(&a, ha, strlen(ha), 0.2, &hit));
  /* a window placed by its offset from the first base of the anchor,
     on either side of it and across it, read the same way whichever
     orientation the read is in */
  const char *reads[2] = {ha, rc};
  for (size_t i = 0; i < 2; i++) {
    const char *r = reads[i];
    size_t len = strlen(r);
    test_assert(anchor_find(&a, r, len, 0.2, &hit));
    test_assert(anchor_bases_window(&hit, r, len, -27, 27, out));
    test_eq_str(out, SITE_NT);
    test_assert(anchor_bases_window(&hit, r, len, 0, 6, out));
    test_eq_str(out, "GGACTG");
    test_assert(anchor_bases_window(&hit, r, len, 24, 6, out));
    test_eq_str(out, "TTTATA");
    /* across the end of the anchor */
    test_assert(anchor_bases_window(&hit, r, len, 21, 6, out));
    test_eq_str(out, "GGCTTT");
    /* past either end of the read */
    test_assert(!anchor_bases_window(&hit, r, len, -40, 6, out));
    test_assert(!anchor_bases_window(&hit, r, len, 30, 12, out));
  }
  /* a read without the anchor */
  const char *other = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
  test_assert(!anchor_find(&a, other, strlen(other), 0.2, &hit));
  /* the bases are read as the read carries them, insertion included */
  const char *inserted = "ATGGAAAAA" "CCTCAAAGAGAGAGAAGAAGAAAAAAGA"
                         FUSION_NT "TTT";
  test_assert(anchor_find(&a, inserted, strlen(inserted), 0.2, &hit));
  test_assert(anchor_bases_before(&hit, inserted, strlen(inserted), 28,
                                  out));
  test_eq_str(out, "CCTCAAAGAGAGAGAAGAAGAAAAAAGA");
  /* the read is searched as it is, as the primers are: a lower-case
     read does not match (the basecallers write upper case) */
  const char *lower = "atggaaaaa" "cctcaaagagagagaagaagaaaaaga"
                      "ggactgtttggggctatagcaggc" "ttt";
  test_assert(!anchor_find(&a, lower, strlen(lower), 0.2, &hit));
  /* no anchor, no read */
  test_assert(!anchor_find(&b, ha, strlen(ha), 0.2, &hit));
  test_assert(!anchor_find(&a, NULL, 0, 0.2, &hit));
  test_assert(anchor_find(&a, ha, strlen(ha), 0.2, &hit));
  test_assert(!anchor_bases_before(&hit, ha, strlen(ha), 0, out));
  test_assert(!anchor_bases_after(&hit, ha, strlen(ha), 0, 0, out));
}
