/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "edlib.h"
#include "alloc.h"
#include "fasta.h"
#include "primer-search.h"
#include "primers.h"
#include "reference.h"
#include "test.h"

struct reference_set_t {
  const fasta_t *fasta;
  size_t probe_len;
  double max_error_rate;
  /* the reverse complement of the probe (probe_len + 1 bytes) */
  char *rc;
};

reference_set_t* reference_set_new(const fasta_t *fasta, size_t probe_len,
    double max_error_rate) {
  assert(fasta != NULL);
  assert(probe_len >= REFERENCE_MIN_PROBE_LEN);
  assert(max_error_rate >= 0 && max_error_rate < 1);
  reference_set_t *self;
  XMALLOC(self, sizeof(reference_set_t));
  self->fasta = fasta;
  self->probe_len = probe_len;
  self->max_error_rate = max_error_rate;
  XMALLOC(self->rc, probe_len + 1);
  return self;
}

void reference_set_delete(reference_set_t *self) {
  if (self == NULL)
    return;
  free(self->rc);
  free(self);
}

size_t reference_set_size(const reference_set_t *self) {
  return self->fasta->n_records;
}

size_t reference_set_probe_len(const reference_set_t *self) {
  return self->probe_len;
}

/* distance of a hit from the terminus from which the read would start */
static long terminus_distance(const reference_hit_t *hit) {
  return hit->reverse ? hit->ref_len - 1 - hit->ref_start : hit->ref_start;
}

/*
   Align <query> (length <len>) to the reference <ref_index> with at most
   <max_errors> errors; if found, and better than <best> (fewer errors,
   or as many and nearer to the terminus, see reference.h), store the
   match in <best>. For <reverse>, <query> is the reverse complement of
   the probe, and the first base of the probe is the last base of the
   alignment.

   (A run of N in a reference matches any probe: the references are
   expected to be complete sequences of the segments.)
*/
static void align_to(const reference_set_t *self, const char *query,
    size_t len, size_t ref_index, bool reverse, long max_errors,
    reference_hit_t *best) {
  const fasta_record_t *ref = self->fasta->records + ref_index;
  if (ref->seq_len < len)
    return;
  EdlibAlignResult result = edlibAlign(query, (int)len, ref->seq,
      (int)ref->seq_len,
      edlibNewAlignConfig((int)max_errors, EDLIB_MODE_HW, EDLIB_TASK_LOC,
                          iupac_equalities, n_iupac_equalities));
  assert(result.status == EDLIB_STATUS_OK);
  if (result.editDistance < 0) {
    edlibFreeAlignResult(result);
    return;
  }
  reference_hit_t hit = {.mapped = true, .ref_index = ref_index,
                         .ref_len = (long)ref->seq_len, .reverse = reverse,
                         .ref_start = -1,
                         .n_errors = result.editDistance};
  /* the location nearest to the terminus (several are reported e.g. for
     a repeat, or for equally good alignments of different lengths) */
  for (int i = 0; i < result.numLocations; i++) {
    long pos = reverse ? result.endLocations[i] : result.startLocations[i];
    if (hit.ref_start < 0 ||
        (reverse ? pos > hit.ref_start : pos < hit.ref_start))
      hit.ref_start = pos;
  }
  edlibFreeAlignResult(result);
  if (!best->mapped || hit.n_errors < best->n_errors ||
      (hit.n_errors == best->n_errors &&
       terminus_distance(&hit) < terminus_distance(best)))
    *best = hit;
}

void reference_set_locate(reference_set_t *self, const char *probe,
    size_t probe_len, reference_hit_t *hit) {
  assert(probe != NULL);
  assert(hit != NULL);
  memset(hit, 0, sizeof(reference_hit_t));
  hit->mapped = false;
  hit->ref_start = -1;
  hit->n_errors = -1;
  size_t len = probe_len < self->probe_len ? probe_len : self->probe_len;
  if (len < REFERENCE_MIN_PROBE_LEN)
    return;
  primers_reverse_complement(probe, len, self->rc);
  self->rc[len] = '\0';
  long max_errors = (long)((double)len * self->max_error_rate);
  for (size_t i = 0; i < self->fasta->n_records; i++) {
    /* a worse match than the best one so far is not needed (an equally
       good one is, for the choice of the location) */
    long k = hit->mapped ? hit->n_errors : max_errors;
    align_to(self, probe, len, i, false, k, hit);
    k = hit->mapped ? hit->n_errors : max_errors;
    align_to(self, self->rc, len, i, true, k, hit);
  }
}

bool reference_hit_is_terminal(const reference_hit_t *hit,
    long terminal_window) {
  assert(hit->mapped);
  return terminus_distance(hit) < terminal_window;
}

/* a deterministic pseudo-random sequence of <len> bases into <seq> */
static void random_seq(char *seq, size_t len, unsigned long *state) {
  for (size_t i = 0; i < len; i++) {
    *state = *state * 6364136223846793005UL + 1442695040888963407UL;
    seq[i] = "ACGT"[(*state >> 33) & 3];
  }
  seq[len] = '\0';
}

/* a reference set of <n> records with the sequences <seqs> */
static fasta_t* test_fasta(char **seqs, size_t n) {
  fasta_t *fa;
  XMALLOC(fa, sizeof(fasta_t));
  XMALLOC(fa->records, sizeof(fasta_record_t) * n);
  fa->n_records = n;
  for (size_t i = 0; i < n; i++) {
    fa->records[i].name = NULL;
    fa->records[i].description = NULL;
    fa->records[i].seq = seqs[i];
    fa->records[i].seq_len = strlen(seqs[i]);
  }
  return fa;
}

/* a probe: <len> bases of <seq> from <pos>, reverse-complemented if
   <rc> (then <pos> is the position of the last base, on the plus strand,
   i.e. the probe is the minus strand read from <pos> backwards) */
static void take_probe(const char *seq, long pos, size_t len, bool rc,
    char *probe) {
  if (rc)
    primers_reverse_complement(seq + pos - (long)len + 1, len, probe);
  else
    memcpy(probe, seq + pos, len);
  probe[len] = '\0';
}

void reference_test(void) {
  enum { REF_LEN = 1000, PROBE = 30 };
  unsigned long state = 42;
  char ref0[REF_LEN + 1], ref1[REF_LEN + 1], probe[64];
  random_seq(ref0, REF_LEN, &state);
  random_seq(ref1, REF_LEN, &state);
  /* a repeat in ref1: the first 40 bases are also at 600 */
  memcpy(ref1 + 600, ref1, 40);
  /* a degenerate base in ref0 */
  ref0[300] = 'R';
  char *seqs[] = {ref0, ref1};
  fasta_t *fa = test_fasta(seqs, 2);
  reference_set_t *rs = reference_set_new(fa, PROBE, 0.2);
  test_eq_size(reference_set_size(rs), 2);
  test_eq_size(reference_set_probe_len(rs), PROBE);
  reference_hit_t hit;

  /* a normal amplicon: the probe is the start of the segment */
  take_probe(ref0, 12, PROBE, false, probe);
  reference_set_locate(rs, probe, PROBE, &hit);
  test_assert(hit.mapped);
  test_eq_size(hit.ref_index, 0);
  test_assert(!hit.reverse);
  test_eq_long(hit.ref_start, 12);
  test_eq_long(hit.n_errors, 0);
  test_eq_long(hit.ref_len, REF_LEN);
  test_assert(reference_hit_is_terminal(&hit, 60));
  test_assert(!reference_hit_is_terminal(&hit, 10));

  /* internal priming: an interior position, with two substitutions */
  take_probe(ref1, 500, PROBE, false, probe);
  probe[5] = probe[5] == 'A' ? 'C' : 'A';
  probe[20] = probe[20] == 'G' ? 'T' : 'G';
  reference_set_locate(rs, probe, PROBE, &hit);
  test_assert(hit.mapped);
  test_eq_size(hit.ref_index, 1);
  test_eq_long(hit.ref_start, 500);
  test_eq_long(hit.n_errors, 2);
  test_assert(!reference_hit_is_terminal(&hit, 60));

  /* an insertion and a deletion in the probe */
  /* (ref0[700..710), an inserted base, ref0[710..720), ref0[721..]) */
  snprintf(probe, sizeof(probe), "%.10s%c%.10s%.*s", ref0 + 700,
           ref0[710] == 'T' ? 'G' : 'T', ref0 + 710, PROBE - 21,
           ref0 + 721);
  test_eq_size(strlen(probe), PROBE);
  reference_set_locate(rs, probe, PROBE, &hit);
  test_assert(hit.mapped);
  test_eq_size(hit.ref_index, 0);
  test_eq_long(hit.ref_start, 700);
  test_assert(hit.n_errors <= 2);

  /* a read of the minus strand from the end of the segment */
  take_probe(ref0, REF_LEN - 1 - 13, PROBE, true, probe);
  reference_set_locate(rs, probe, PROBE, &hit);
  test_assert(hit.mapped);
  test_assert(hit.reverse);
  test_eq_long(hit.ref_start, REF_LEN - 1 - 13);
  test_eq_long(hit.n_errors, 0);
  test_assert(reference_hit_is_terminal(&hit, 60));

  /* a read of the minus strand from an interior position */
  take_probe(ref1, 429, PROBE, true, probe);
  reference_set_locate(rs, probe, PROBE, &hit);
  test_assert(hit.mapped);
  test_assert(hit.reverse);
  test_eq_size(hit.ref_index, 1);
  test_eq_long(hit.ref_start, 429);
  test_assert(!reference_hit_is_terminal(&hit, 60));

  /* a repeat: equally good at the start and inside, the terminal
     location is chosen (an ambiguous read is kept intact) */
  take_probe(ref1, 5, PROBE, false, probe);
  reference_set_locate(rs, probe, PROBE, &hit);
  test_assert(hit.mapped);
  test_eq_long(hit.ref_start, 5);
  test_assert(reference_hit_is_terminal(&hit, 60));

  /* a degenerate base in the reference matches its bases */
  take_probe(ref0, 290, PROBE, false, probe);
  probe[10] = 'G';
  reference_set_locate(rs, probe, PROBE, &hit);
  test_eq_long(hit.ref_start, 290);
  test_eq_long(hit.n_errors, 0);

  /* a longer probe is cut to the probe length of the set: the bases
     after it are not used */
  take_probe(ref1, 250, PROBE + 20, false, probe);
  memset(probe + PROBE, 'A', 20);
  reference_set_locate(rs, probe, PROBE + 20, &hit);
  test_eq_long(hit.ref_start, 250);
  test_eq_long(hit.n_errors, 0);

  /* a sequence not in the references */
  random_seq(probe, PROBE, &state);
  reference_set_locate(rs, probe, PROBE, &hit);
  test_assert(!hit.mapped);
  test_eq_long(hit.ref_start, -1);

  /* a probe too short to be located */
  take_probe(ref1, 500, REFERENCE_MIN_PROBE_LEN - 1, false, probe);
  reference_set_locate(rs, probe, REFERENCE_MIN_PROBE_LEN - 1, &hit);
  test_assert(!hit.mapped);

  reference_set_delete(rs);
  free(fa->records);
  free(fa);
}
