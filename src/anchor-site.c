/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "anchor-site.h"
#include "codons.h"
#include "edlib.h"
#include "iupac.h"
#include "test.h"

static char complement(char c) {
  switch (c) {
    case 'A': case 'a': return 'T';
    case 'C': case 'c': return 'G';
    case 'G': case 'g': return 'C';
    case 'T': case 't': case 'U': case 'u': return 'A';
    default: return 'N';
  }
}

static char upper(char c) {
  return c >= 'a' && c <= 'z' ? (char)(c - 'a' + 'A') : c;
}

/* write into <out> the reverse complement of the <len> characters of
   <seq>, terminated */
static void revcomp(const char *seq, size_t len, char *out) {
  for (size_t i = 0; i < len; i++)
    out[i] = complement(seq[len - 1 - i]);
  out[len] = '\0';
}

bool anchor_of_peptide(anchor_t *self, size_t ref, const char *seq,
    size_t seq_len, const char *peptide, size_t len) {
  assert(self != NULL);
  memset(self, 0, sizeof(anchor_t));
  if (len == 0 || len > ANCHOR_MAX_LEN)
    return false;
  long pos = codons_find_peptide(seq, seq_len, peptide);
  if (pos < 0 || (size_t)pos + len > seq_len)
    return false;
  self->found = true;
  self->ref = ref;
  self->ref_pos = pos;
  self->len = len;
  for (size_t i = 0; i < len; i++)
    self->seq[i] = upper(seq[pos + i]);
  self->seq[len] = '\0';
  revcomp(self->seq, len, self->rev);
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
  EdlibAlignResult r = edlibAlign(query, (int)query_len, target,
      (int)target_len,
      edlibNewAlignConfig((int)max_errors, EDLIB_MODE_HW, EDLIB_TASK_LOC,
                          iupac_equalities, n_iupac_equalities));
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

bool anchor_bases_before(const anchor_t *self, const char *read,
    size_t read_len, double max_error_rate, size_t n, char *out) {
  assert(self != NULL);
  assert(out != NULL);
  if (!self->found || n == 0 || read == NULL || read_len == 0)
    return false;
  long max_errors = (long)(max_error_rate * (double)self->len);
  long fwd_start = 0, fwd_end = 0, rev_start = 0, rev_end = 0;
  long fwd = search(self->seq, self->len, read, read_len, max_errors,
                    &fwd_start, &fwd_end);
  long rev = search(self->rev, self->len, read, read_len, max_errors,
                    &rev_start, &rev_end);
  if (fwd < 0 && rev < 0)
    return false;
  /* as good in both orientations: which side of the match the bases are
     on is not decided */
  if (fwd >= 0 && rev >= 0 && fwd == rev)
    return false;
  if (rev < 0 || (fwd >= 0 && fwd < rev)) {
    /* the read runs as the reference: the bases are those before the
       match */
    if (fwd_start < (long)n)
      return false;
    for (size_t i = 0; i < n; i++)
      out[i] = upper(read[fwd_start - (long)n + (long)i]);
    out[n] = '\0';
    return true;
  }
  /* the read is the reverse complement of the reference: the bases
     before the anchor there are those after the match here, reverse
     complemented */
  if (rev_end + (long)n >= (long)read_len)
    return false;
  revcomp(read + rev_end + 1, n, out);
  return true;
}

/* --- unit tests --- */

/* the fusion peptide GLFGAIAG, which begins HA2 */
#define FUSION_NT "GGACTGTTTGGGGCTATAGCAGGC"
/* the cleavage site before it in an H5 of clade 2.3.4.4b: PQRERRRKR */
#define SITE_NT "CCTCAAAGAGAGAGAAGAAGAAAAAGA"

void anchor_site_test(void) {
  const char *ha = "ATGGAAAAA" SITE_NT FUSION_NT "TTTATAGAGGGA";
  anchor_t a;
  test_assert(anchor_of_peptide(&a, 3, ha, strlen(ha), "GLFGA", 24));
  test_assert(a.found);
  test_eq_size(a.ref, 3);
  test_eq_long(a.ref_pos, 9 + (long)strlen(SITE_NT));
  test_eq_str(a.seq, FUSION_NT);
  test_eq_str(a.rev, "GCCTGCTATAGCCCCAAACAGTCC");
  /* a peptide which is not there, and an anchor past the end */
  anchor_t b;
  test_assert(!anchor_of_peptide(&b, 0, ha, strlen(ha), "WWWWW", 24));
  test_assert(!anchor_of_peptide(&b, 0, ha, strlen(ha), "GLFGA", 48));
  test_assert(!anchor_of_peptide(&b, 0, ha, strlen(ha), "GLFGA", 0));

  char out[32];
  /* the read as the reference: the last 9 codons before HA2 */
  test_assert(anchor_bases_before(&a, ha, strlen(ha), 0.2, 27, out));
  test_eq_str(out, SITE_NT);
  /* the same read, reverse complemented */
  char rc[128];
  revcomp(ha, strlen(ha), rc);
  test_assert(anchor_bases_before(&a, rc, strlen(rc), 0.2, 27, out));
  test_eq_str(out, SITE_NT);
  /* a read with errors in the anchor still finds it */
  char noisy[128];
  snprintf(noisy, sizeof(noisy), "ATGGAAAAA" SITE_NT
           "GGACTCTTTGGGGCTTTAGCAGGC" "TTTATAGAGGGA");
  test_assert(anchor_bases_before(&a, noisy, strlen(noisy), 0.2, 27, out));
  test_eq_str(out, SITE_NT);
  /* a read which does not reach back far enough */
  const char *short_read = "AAAGA" FUSION_NT "TTT";
  test_assert(!anchor_bases_before(&a, short_read, strlen(short_read), 0.2,
                                   27, out));
  test_assert(anchor_bases_before(&a, short_read, strlen(short_read), 0.2,
                                  5, out));
  test_eq_str(out, "AAAGA");
  /* a read without the anchor */
  const char *other = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
  test_assert(!anchor_bases_before(&a, other, strlen(other), 0.2, 27, out));
  /* the bases are read as the read carries them, insertion included */
  const char *inserted = "ATGGAAAAA" "CCTCAAAGAGAGAGAAGAAGAAAAAAGA"
                         FUSION_NT "TTT";
  test_assert(anchor_bases_before(&a, inserted, strlen(inserted), 0.2, 28,
                                  out));
  test_eq_str(out, "CCTCAAAGAGAGAGAAGAAGAAAAAAGA");
  /* the read is searched as it is, as the primers are: a lower-case
     read does not match (the basecallers write upper case) */
  const char *lower = "atggaaaaa" "cctcaaagagagagaagaagaaaaaga"
                      "ggactgtttggggctatagcaggc" "ttt";
  test_assert(!anchor_bases_before(&a, lower, strlen(lower), 0.2, 27, out));
  /* no anchor, no read */
  test_assert(!anchor_bases_before(&b, ha, strlen(ha), 0.2, 27, out));
  test_assert(!anchor_bases_before(&a, NULL, 0, 0.2, 27, out));
  test_assert(!anchor_bases_before(&a, ha, strlen(ha), 0.2, 0, out));
}
