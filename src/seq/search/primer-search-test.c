/* (c) 2026 Giorgio Gonnella */

/*
   Unit tests of the search of a primer in a read (primer-search.h): in
   the whole read, and at an end of it with outer bases missing.
*/

#include <stdbool.h>
#include <stddef.h>
#include <string.h>
#include "primer-search.h"
#include "primer-search-test.h"
#include "primers.h"
#include "test.h"

/* run a search and compare the result with the expected values
   (start == -1: not found) */
static void check(primer_search_t *s, const char *seq, const char *read,
    long max_errors, bool left, bool right, long start, long end,
    bool partial, long n_errors) {
  primer_hit_t hit;
  bool found = primer_search_run_seq(s, seq, strlen(seq), read,
      strlen(read), max_errors, left, right, &hit);
  test_assert(found == (start >= 0));
  test_eq(hit.start, start, long, "%ld");
  test_eq(hit.end, end, long, "%ld");
  if (found) {
    test_assert(hit.partial == partial);
    test_eq(hit.n_errors, n_errors, long, "%ld");
  }
}

/* the bounds of the terminus of a primer found in a read */
static void check_terminus(primer_search_t *s, const char *name,
    const char *read, long max_errors, bool left, bool right,
    long terminus_start, long terminus_end) {
  const primer_t *p = primers_builtin(name);
  test_assert(p != NULL);
  primer_hit_t hit;
  test_assert(primer_search_run(s, p, read, strlen(read), max_errors, left,
                                right, &hit));
  test_eq(hit.terminus_start, terminus_start, long, "%ld");
  test_eq(hit.terminus_end, terminus_end, long, "%ld");
}

/* the costs of primer_search_end() for a primer of length <len>: up to
   <max_skip> outer bases missing at the cost <skip_c>, up to <max_edge>
   at the end of the read at the cost <edge_c> (or <skip_c>, if lower) */
static void end_costs(size_t len, size_t max_skip, long skip_c,
    size_t max_edge, long edge_c, long *skip, long *edge) {
  for (size_t q = 0; q <= len; q++) {
    skip[q] = q == 0 ? 0 : (q <= max_skip ? skip_c : -1);
    edge[q] = q == 0 ? 0 : (q <= max_edge ? edge_c : -1);
    if (skip[q] >= 0 && (edge[q] < 0 || skip[q] < edge[q]))
      edge[q] = skip[q];
  }
}

/* run a search at an end and compare the result with the expected
   values (<start> == -1: not found) */
static void check_end(primer_search_t *s, const primer_t *p,
    const char *read, bool left, const long *skip, const long *edge,
    long max_cost, long start, long end, long n_errors,
    long terminus_start, long terminus_end) {
  primer_hit_t hit;
  bool found = primer_search_end(s, p, read, strlen(read), left, skip,
                                 edge, max_cost, &hit);
  test_assert(found == (start >= 0));
  test_eq(hit.start, start, long, "%ld");
  test_eq(hit.end, end, long, "%ld");
  if (found) {
    test_eq(hit.n_errors, n_errors, long, "%ld");
    test_eq(hit.terminus_start, terminus_start, long, "%ld");
    test_eq(hit.terminus_end, terminus_end, long, "%ld");
  }
}

/* the search at an end of the read, with missing outer primer bases */
static void test_end(primer_search_t *s) {
  long skip[32], edge[32];
  char read[256];
  /* mbtuni-12: 10 bases of tail, then the 12 bases of the terminus */
  const primer_t *p12 = primers_builtin("mbtuni-12");
  const char *body = "GGTTCCAAGGTTCCAAGGTTCCAA";
  /* up to 11 outer bases missing, at the cost of one error */
  end_costs(p12->seq_len, 11, 1, 11, 1, skip, edge);
  /* the whole primer after a flank */
  snprintf(read, sizeof(read), "TTTTT%s%s", p12->seq, body);
  check_end(s, p12, read, true, skip, edge, 2, 5, 26, 0, 15, 27);
  /* without its first 3 bases: the errors do not include the cost */
  snprintf(read, sizeof(read), "TTTTT%s%s", p12->seq + 3, body);
  check_end(s, p12, read, true, skip, edge, 2, 5, 23, 0, 12, 24);
  /* without 11 bases, one of the terminus: the terminus begins with
     the alignment */
  snprintf(read, sizeof(read), "TTTTT%s%s", p12->seq + 11, body);
  check_end(s, p12, read, true, skip, edge, 2, 5, 15, 0, 5, 16);
  /* 12 bases missing: not allowed (with at most one error, as the
     flank base before them is not the missing primer base) */
  snprintf(read, sizeof(read), "TTTTT%s%s", p12->seq + 12, body);
  check_end(s, p12, read, true, skip, edge, 1, -1, -1, 0, 0, 0);
  /* the read begins inside the primer: with a cost of 0 at the end of
     the read, the primer overhangs it, and the terminus begins one base
     before the read */
  end_costs(p12->seq_len, 11, 1, 11, 0, skip, edge);
  snprintf(read, sizeof(read), "%s%s", p12->seq + 11, body);
  check_end(s, p12, read, true, skip, edge, 0, 0, 10, 0, -1, 11);
  /* a read with an insertion and one with a deletion in the tail (as
     with primer_search_run(), the terminus moves with them) */
  end_costs(p12->seq_len, 0, 0, 0, 0, skip, edge);
  check_end(s, p12, "ACGCGTGATCCAGCAAAAGCAGGTTTTTTTTT", true, skip, edge,
            2, 0, 22, 1, 11, 23);
  check_end(s, p12, "ACGCGTGAT" "AGCAAAAGCAGGTTTTTTTTTT", true, skip, edge,
            2, 0, 20, 1, 9, 21);

  /*
     The cost of the missing bases: a shorter match is allowed fewer
     errors. Uni12 without its first two bases (cost 1) after a flank,
     with one mismatch, is found; with two, it is not (the whole primer
     would be allowed two errors).
  */
  const primer_t *u12 = primers_builtin("uni12");
  end_costs(u12->seq_len, 4, 1, 4, 1, skip, edge);
  check_end(s, u12, "CCCCCCCAAATGCAGGTTCCAAGGTTCCAA", true, skip, edge, 2,
            6, 15, 1, 6, 16);
  check_end(s, u12, "CCCCCCCAAATGCTGGTTCCAAGGTTCCAA", true, skip, edge, 2,
            -1, -1, 0, 0, 0);
  /* the lowest cost is taken: the inner exact primer rather than an
     outer one with an error */
  snprintf(read, sizeof(read), "TTAGCAAATGCAGG%sAGCAAAAGCAGG%s", body,
           body);
  check_end(s, u12, read, true, skip, edge, 2, 38, 49, 0, 38, 50);
  /* ... and among equal costs the outermost: the truncated primer (cost
     1, no error) rather than an inner match with one error */
  snprintf(read, sizeof(read), "TTCAAAAGCAGG%sAGCAAATGCAGG%s", body, body);
  check_end(s, u12, read, true, skip, edge, 2, 2, 11, 0, 2, 12);
  /*
     Ends of equal cost: the outermost, so that a base which may be
     template is kept. A simulated read whose Uni12 has two errors, the
     second at its last base (G read as A), followed by the template G
     of segment position 13: the alignments ending at the G before the
     A (the last primer base deleted), at the A and at the template G
     (the A inserted) all cost two errors. The template G is kept, and
     so is the A.
  */
  end_costs(u12->seq_len, 0, 0, 0, 0, skip, edge);
  check_end(s, u12, "TTTTTAGCAAATGCAGAGTAGATATTCCAAGG", true, skip, edge, 2,
            5, 15, 2, 5, 16);
  /* the price: the last primer base, with an error (G read as T), is
     kept when it can also be a template base after the G deleted */
  check_end(s, u12, "TTTTTAGCAAAAGCAGTCCTTAAGGTTCCAA", true, skip, edge, 2,
            5, 15, 1, 5, 16);
  /* degenerate bases: the R of the primer, an N of the read */
  const primer_t *u12r = primers_builtin("uni12r");
  end_costs(u12r->seq_len, 0, 0, 0, 0, skip, edge);
  check_end(s, u12r, "TTAGCGAAAGCAGGTTCCAAGGTTCCAA", true, skip, edge, 0,
            2, 13, 0, 2, 14);
  check_end(s, u12r, "TTAGCAAANGCAGGTTCCAAGGTTCCAA", true, skip, edge, 0,
            2, 13, 0, 2, 14);

  /*
     The right end: the reverse complement of mbtuni-13, whose terminus
     (13 bases) is its prefix; the missing bases are its last ones
  */
  const primer_t *p13 = primers_builtin("mbtuni-13");
  // (0-terminated, as it is printed)
  char rc_seq[32] = "";
  primer_t rc;
  primers_reverse(p13, rc_seq, &rc);
  long n = (long)strlen(body);
  end_costs(rc.seq_len, 11, 1, 11, 0, skip, edge);
  snprintf(read, sizeof(read), "%s%sTTTTT", body, rc.seq);
  check_end(s, &rc, read, false, skip, edge, 2, n, n + 22, 0, n, n + 13);
  /* without its last 3 bases, before a flank */
  snprintf(read, sizeof(read), "%s%.20sTTTTT", body, rc.seq);
  check_end(s, &rc, read, false, skip, edge, 2, n, n + 19, 0, n, n + 13);
  /* without 11 bases, one of the terminus: it ends with the alignment */
  snprintf(read, sizeof(read), "%s%.12sGGGGG", body, rc.seq);
  check_end(s, &rc, read, false, skip, edge, 2, n, n + 11, 0, n, n + 12);
  /* the read ends inside the terminus: it ends one base after the read */
  snprintf(read, sizeof(read), "%s%.12s", body, rc.seq);
  check_end(s, &rc, read, false, skip, edge, 0, n, n + 11, 0, n, n + 13);

  /*
     A real Nanopore read of the Opti1 primers (SRR37040770) primed at
     the internal site of HA nucleotide 1087 at both ends: its right end
     is the reverse complement of opti1-f2 (tail and the 12 bases of
     terminus-12-g) after the template base C (G on the plus strand).
     The right primer of the scheme, opti1-r1, matches it with four
     errors (its terminus has 13 bases) both with and without that C:
     the C is kept.
  */
  const primer_t *r1 = primers_builtin("opti1-r1");
  char rc1_seq[32] = "";
  primer_t rc1;
  primers_reverse(r1, rc1_seq, &rc1);
  end_costs(rc1.seq_len, 0, 0, 0, 0, skip, edge);
  snprintf(read, sizeof(read), "%sC%s%sTTTTT", body, "CCTGCTTTCGCT",
           "GGCGCGTAAC");
  check_end(s, &rc1, read, false, skip, edge, 4, n + 1, n + 22, 4, n + 1,
            n + 13);
}

void primer_search_test(void) {
  primer_search_t *s = primer_search_new();
  /* sequences without degenerate bases */
  check(s, "AGCTTA", "GGAGCTTAGG", 1, false, false, 2, 7, false, 0);
  check(s, "AGCTTA", "GGAGGTTAGG", 1, false, false, 2, 7, false, 1);
  check(s, "ACGTACGT", "", 2, false, false, -1, -1, false, 0);
  check(s, "ACGTACGT", "TACGT", 2, true, false, 0, 4, true, 0);
  check(s, "ACGTACGT", "GGACGTAC", 2, false, true, 2, 7, false, 2);
  check(s, "ACGTACGT", "GGGGGGACGTA", 2, false, true, 6, 10, true, 0);
  check(s, "ACGTACGT", "GGGGGGACGTA", 2, false, false, -1, -1, false, 0);
  /* degenerate bases of the primer (R = A or G) match either base */
  check(s, "ACRTACGT", "GGACATACGTGG", 0, false, false, 2, 9, false, 0);
  check(s, "ACRTACGT", "GGACGTACGTGG", 0, false, false, 2, 9, false, 0);
  check(s, "ACRTACGT", "GGACCTACGTGG", 0, false, false, -1, -1, false, 0);
  /* a degenerate base of the read (N) matches any base of the primer */
  check(s, "ACGTACGT", "GGACNTACGTGG", 0, false, false, 2, 9, false, 0);
  /* ... and any degenerate base of the primer, also when the primer
     overhangs the read (a read starting inside the degenerate Uni12
     primer AGCRAAAGCAGG, after its R, is found without errors) */
  check(s, "ACRTACGT", "GGACNTACGTGG", 0, false, false, 2, 9, false, 0);
  check(s, "AGCRAAAGCAGG", "AAGCAGGTTTTTTTTT", 0, true, false, 0, 6, true,
        0);

  /*
     The terminus of a tailed primer: mbtuni-12 is
     ACGCGTGATC + AGCAAAAGCAGG (10 bases of tail, then the 12 bases of
     the conserved terminus, see primers.h)
  */
  /* the whole primer at the start of the read: the terminus begins at 10 */
  check_terminus(s, "mbtuni-12", "ACGCGTGATCAGCAAAAGCAGGTTTTTTTTTT", 2,
                 true, false, 10, 22);
  /* the primer after two bases of leader */
  check_terminus(s, "mbtuni-12", "TTACGCGTGATCAGCAAAAGCAGGTTTTTTTT", 2,
                 true, false, 12, 24);
  /* a deletion in the tail moves the terminus one base to the left */
  check_terminus(s, "mbtuni-12", "ACGCGTGAT" "AGCAAAAGCAGGTTTTTTTTTT", 2,
                 true, false, 9, 21);
  /* an insertion in the tail moves it one base to the right */
  check_terminus(s, "mbtuni-12", "ACGCGTGATCCAGCAAAAGCAGGTTTTTTTTT", 2,
                 true, false, 11, 23);
  /* the read begins inside the terminus (the primer overhangs the read
     on the left, here by half of its length, the maximum): the terminus
     begins one base before the read */
  check_terminus(s, "mbtuni-12", "GCAAAAGCAGGTTTTTTTTTTTTTTTTTTTT", 2,
                 true, false, -1, 11);
  /* a longer overhang (only 9 of the 22 bases of the primer are in the
     read) is not accepted */
  {
    const primer_t *p = primers_builtin("mbtuni-12");
    primer_hit_t hit;
    test_assert(!primer_search_run(s, p, "AAAAGCAGGTTTTTTTTTTTTTTTTTTTT",
                                   29, 2, true, false, &hit));
  }
  /* a primer without a tail: the terminus is the whole alignment */
  check_terminus(s, "uni12", "TTAGCAAAAGCAGGTTTTTTTTTTTTTT", 2, true,
                 false, 2, 14);
  /* the reverse complement of a primer, as it is found at the other end
     of an amplicon: its terminus is its prefix */
  {
    const primer_t *p = primers_builtin("mbtuni-13");
    char rc_seq[32];
    primer_t rc;
    primers_reverse(p, rc_seq, &rc);
    test_eq_size(rc.terminus_from, 0);
    test_eq_mem(rc.seq, rc.seq_len, "CCTTGTTTCTACTGATCACGCGT");
    primer_hit_t hit;
    /* 10 bases of segment, the terminus, then the reverse tail */
    const char *read = "TTTTTTTTTT" "CCTTGTTTCTACT" "GATCACGCGT";
    test_assert(primer_search_run(s, &rc, read, strlen(read), 2, false,
                                  true, &hit));
    test_eq(hit.start, 10, long, "%ld");
    test_eq(hit.terminus_start, 10, long, "%ld");
    test_eq(hit.terminus_end, 23, long, "%ld");
  }
  test_end(s);
  primer_search_delete(s);
}
