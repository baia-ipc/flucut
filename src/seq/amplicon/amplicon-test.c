/* (c) 2026 Giorgio Gonnella */

/*
   Unit tests of the analysis of a read (amplicon.h): the primers at the
   ends, the trimming policies, internal priming, the adapters, the
   classification against the references and the ends located from
   inside.
*/

#include <stdio.h>
#include <string.h>
#include "adapters.h"
#include "amplicon.h"
#include "amplicon-ends.h"
#include "amplicon-test.h"
#include "fasta.h"
#include "logger.h"
#include "reference.h"
#include "scheme.h"
#include "sequtil.h"
#include "test.h"

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

/* a deterministic pseudo-random sequence of <len> bases into <seq> */
static void test_random_seq(char *seq, size_t len, unsigned long state) {
  for (size_t i = 0; i < len; i++) {
    state = state * 6364136223846793005UL + 1442695040888963407UL;
    seq[i] = "ACGT"[(state >> 33) & 3];
  }
  seq[len] = '\0';
}

/*
   The ends located from inside the read: a sample whose first bases
   after the Uni12 primer differ from the reference, so that the probe
   of the left end does not map
*/
static void test_inner_ends(Logger *logger) {
  enum { INTERIOR = 900 };
  const char *uni12 = "AGCAAAAGCAGG";
  const char *rc13 = "CCTTGTTTCTACTGATCACGCGT";       /* rc mbtuni-13 */
  char interior[INTERIOR + 1], sample[INTERIOR + 1], segment[1024];
  test_random_seq(interior, INTERIOR, 7);
  memcpy(sample, interior, INTERIOR + 1);
  for (size_t i = 0; i < 30; i += 3)
    sample[i] = sample[i] == 'A' ? 'C' : 'A';
  snprintf(segment, sizeof(segment), "%s%sCCTTGTTTCTACT", uni12, interior);
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
  char read[1024];

  /* a complete amplicon: the left end is not located by its probe; it
     is by a probe further inside, at the terminus; the origin of the
     end, and so the trimming, stay those of an end not located */
  snprintf(read, sizeof(read), "%s%s%s", uni12, sample, rc13);
  const amplicon_t *r = amplicon_analyzer_run(a, read, strlen(read));
  test_assert(r->left.found && !r->left.ref.mapped);
  test_assert(r->left.origin == ORIGIN_UNKNOWN);
  test_eq_long(r->keep_from, 12);
  test_assert(r->left.inner.mapped);
  test_eq_long(r->left.inner.ref_start, 12);
  test_assert(amplicon_end_position(&r->right) == &r->right.ref);
  test_assert(!r->junction.found);

  /* an internal deletion with a left flank of 50 bases: the probe which
     is located is of the right flank, beyond the junction; the bases
     skipped do not continue the reference there, and agree with it at
     the terminus: the end is at the terminus, and the junction is found */
  snprintf(read, sizeof(read), "%s%.50s%s%s", uni12, sample, sample + 600,
           rc13);
  r = amplicon_analyzer_run(a, read, strlen(read));
  test_assert(!r->left.ref.mapped && r->left.inner.mapped);
  test_eq_long(r->left.inner.ref_start, 12);
  test_assert(r->junction.found);
  test_assert(r->junction.left_to >= 60 && r->junction.left_to <= 64);
  test_eq_long(r->junction.right_from - r->junction.left_to, 550);

  /* a fragment without primers, with an internal deletion: both ends
     are located from inside */
  snprintf(read, sizeof(read), "%.200s%.200s", sample + 100, sample + 700);
  r = amplicon_analyzer_run(a, read, strlen(read));
  test_eq_size(amplicon_n_ends(r), 0);
  test_assert(r->left.inner.mapped && r->right.inner.mapped);
  test_eq_long(r->left.inner.ref_start, 112);
  test_eq_long(r->right.inner.ref_start, 911);
  test_assert(r->junction.found);
  /* (bases which fit both flanks are given to the right one) */
  test_assert(r->junction.left_to >= 308 && r->junction.left_to <= 312);
  test_eq_long(r->junction.right_from - r->junction.left_to, 400);

  amplicon_analyzer_delete(a);
  scheme_delete(sc);
  reference_set_delete(refs);
  reference_index_delete(index);
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
  amplicon_missing_costs(0.2, 12, 2, skip, edge);
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
  /* the left primer at both ends, at the end of the read
     reverse-complemented (two amplicons joined without a primer between
     them): both are removed (with the untailed Uni12 primer: the tailed
     ones share their tail, see amplicon_analyze()) */
  scheme_t *untailed = scheme_new_builtin(SCHEME_DEFAULT, logger);
  amplicon_analyzer_t *u = amplicon_analyzer_new(untailed, 0.2,
      TRIM_PRIMER, AMPLICON_DEFAULT_WINDOW);
  snprintf(read, sizeof(read), "AGCAAAAGCAGG%sCCTGCTTTTGCT", body);
  check_keep(u, read, STRAND_PLUS, 1, body);
  amplicon_analyzer_delete(u);
  scheme_delete(untailed);
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
  test_inner_ends(logger);
  delete_logger(logger);
}
