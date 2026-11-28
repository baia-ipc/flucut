/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cleavage-counts-test.h"
#include "cleavage-counts.h"
#include "codons.h"
#include "examples-table.h"
#include "test.h"

/* --- unit tests --- */

/* the cleavage site of a highly pathogenic and of a low pathogenic
   haemagglutinin, as the reads carry them (metadata/examples.yaml) */
#define HPAI_NT EXAMPLE_CLEAVAGE_SITE_MULTIBASIC
#define LPAI_NT EXAMPLE_CLEAVAGE_SITE_MONOBASIC

static void add_n(cleavage_counts_t *c, const char *nt, size_t n) {
  for (size_t i = 0; i < n; i++) {
    cleavage_counts_searched(c);
    cleavage_counts_add(c, nt);
  }
}

static void test_counts(void) {
  cleavage_counts_t c;
  cleavage_counts_clear(&c);
  test_eq_size(cleavage_counts_reads(&c), 0);
  cleavage_result_t r;
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_assert(!r.found);
  test_eq_size(r.verdict, CLEAVAGE_TOO_FEW_READS);

  /* 70 reads of the motif, 20 of a spelling of the same residues (the
     last codon AGG instead of AGA, both R) and 10 whose last residue is
     an isoleucine instead */
  add_n(&c, HPAI_NT, 70);
  add_n(&c, "CAAAGAGAGAGAAGAAGAAAAAGG", 20);
  add_n(&c, "CAAAGAGAGAGAAGAAGAAAAATA", 10);
  test_eq_size(cleavage_counts_reads(&c), 100);
  test_eq_size(c.n_motifs, 3);
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_assert(r.found);
  test_eq_str(r.aa, "QRERRRKR");
  /* the most frequent spelling of those residues */
  test_eq_str(r.nt, HPAI_NT);
  test_eq_size(r.reads, 100);
  /* the two spellings of the same residues are counted together */
  test_eq_size(r.motif_reads, 90);
  test_eq_size(r.distinct, 3);
  test_eq_size(r.distinct_proteins, 2);
  /* both readings are carried by at least a tenth of the reads */
  test_eq_size(r.supported_proteins, 2);
  test_eq_size(r.basic, 5);
  test_eq_size(r.reading, CLEAVAGE_MULTIBASIC);
  test_eq_size(r.verdict, CLEAVAGE_FIRM);
  test_assert(r.support > 0.89 && r.support < 0.91);
  test_assert(r.multibasic_share > 0.89 && r.multibasic_share < 0.91);
  /* without a reference site nothing is compared with it */
  test_eq_size(r.reference_reading, CLEAVAGE_N_READINGS);
  test_assert(!r.disagrees_with_reference);
  /* the sites found, most frequent first */
  test_eq_size(r.n_top, 2);
  test_eq_str(r.top[0].aa, "QRERRRKR");
  test_eq_size(r.top[0].reads, 90);
  test_eq_size(r.top[0].reading, CLEAVAGE_MULTIBASIC);
  test_eq_str(r.top[1].aa, "QRERRRKI");
  test_eq_size(r.top[1].reads, 10);
  test_eq_size(r.top[1].reading, CLEAVAGE_UNREADABLE);

  /* the reads read like the reference: still firm */
  cleavage_counts_set_reference(&c, "LREKRRKR", "GLFGAIAG");
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_str(r.reference_aa, "LREKRRKR");
  test_eq_size(r.reference_reading, CLEAVAGE_MULTIBASIC);
  test_assert(!r.disagrees_with_reference);
  test_eq_size(r.verdict, CLEAVAGE_FIRM);
  /* a sample whose reference is multibasic but whose reads read as a
     low pathogenic site holds the reads of another haemagglutinin: it
     is never firm, however many reads agree */
  cleavage_counts_t mixed;
  cleavage_counts_clear(&mixed);
  cleavage_counts_set_reference(&mixed, "LREKRRKR", "GLFGAIAG");
  add_n(&mixed, LPAI_NT, 100);
  cleavage_result_of(&mixed, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_str(r.aa, "NVPARSSR");
  test_assert(r.support > 0.99);
  test_assert(r.disagrees_with_reference);
  test_eq_size(r.verdict, CLEAVAGE_CONFIRM);

  /* the verdict does not depend on the chemistry: the same reads are
     judged the same whatever they were sequenced with */
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R9, &r);
  test_eq_size(r.verdict, CLEAVAGE_FIRM);
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_ILLUMINA, &r);
  test_eq_size(r.verdict, CLEAVAGE_FIRM);

  /* too few reads */
  cleavage_counts_t few;
  cleavage_counts_clear(&few);
  add_n(&few, HPAI_NT, CLEAVAGE_MIN_READS - 1);
  cleavage_result_of(&few, CLEAVAGE_CHEMISTRY_R10, &r);
  test_assert(r.found);
  test_eq_size(r.verdict, CLEAVAGE_TOO_FEW_READS);
  add_n(&few, HPAI_NT, 1);
  cleavage_result_of(&few, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_size(r.verdict, CLEAVAGE_FIRM);

  /* merging adds the counts of both */
  cleavage_counts_t other;
  cleavage_counts_clear(&other);
  add_n(&other, HPAI_NT, 5);
  add_n(&other, LPAI_NT, 1);
  cleavage_counts_merge(&c, &other);
  test_eq_size(cleavage_counts_reads(&c), 106);
  test_eq_size(c.searched, 106);
  test_eq_size(c.n_motifs, 4);
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_size(r.motif_reads, 95);
}

/* a reading of the site carried by many reads and by a share of them,
   and readings carried by one read each: the raw counts include the
   single reads, the count of supported readings does not */
static void test_supported_readings(void) {
  cleavage_counts_t c;
  cleavage_counts_clear(&c);
  add_n(&c, HPAI_NT, 60);
  add_n(&c, LPAI_NT, 30);
  /* ten reads with a base of their own changed in the first 10 codon
     positions after the first: ten spellings, one read each */
  char err[CLEAVAGE_SITE_NT + 1];
  for (size_t i = 0; i < 10; i++) {
    memcpy(err, HPAI_NT, sizeof(err));
    err[3 + i] = err[3 + i] == 'T' ? 'C' : 'T';
    add_n(&c, err, 1);
  }
  cleavage_result_t r;
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_size(r.reads, 100);
  test_eq_size(r.distinct, 12);
  /* the readings of the single reads are counted here ... */
  test_assert(r.distinct_proteins > 3);
  /* ... and not here: the two readings of the sample */
  test_eq_size(r.supported_proteins, 2);
  test_eq_str(r.top[0].aa, "QRERRRKR");
  test_eq_size(r.top[1].reads, 30);

  /* a reading carried by two reads, below the share, is not supported;
     at the share it is */
  cleavage_counts_t d;
  cleavage_counts_clear(&d);
  add_n(&d, HPAI_NT, 98);
  add_n(&d, LPAI_NT, 2);
  cleavage_result_of(&d, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_size(r.supported_proteins, 1);
  add_n(&d, LPAI_NT, 10);
  cleavage_result_of(&d, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_size(r.supported_proteins, 2);
}

/* the control windows: reads which agree at the site but nowhere else
   are not firm */
static void test_control(void) {
  cleavage_counts_t c;
  cleavage_counts_clear(&c);
  add_n(&c, HPAI_NT, 100);
  cleavage_result_t r;
  /* no control window read: nothing is known about the rest, and the
     site alone decides */
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_size(r.control_reads, 0);
  test_assert(!r.control_weak);
  test_eq_size(r.verdict, CLEAVAGE_FIRM);
  /* control windows the reads agree on (the second motif differs in a
     residue, not only in a base: the control counts residues too) */
  for (size_t i = 0; i < 90; i++)
    for (size_t w = 0; w < CLEAVAGE_N_CONTROLS; w++)
      cleavage_counts_add_control(&c, w, "TTTATAGAGGGAGGATGGCAGGGA");
  for (size_t i = 0; i < 10; i++)
    for (size_t w = 0; w < CLEAVAGE_N_CONTROLS; w++)
      cleavage_counts_add_control(&c, w, "TTTATAGAGGGAGGATGGCAGTCA");
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_size(r.control_reads, 200);
  test_assert(r.control_support > 0.89 && r.control_support < 0.91);
  test_assert(!r.control_weak);
  test_eq_size(r.verdict, CLEAVAGE_FIRM);
  /* reads which agree nowhere, not even where nothing is difficult:
     the site is not wrong so much as unknown */
  cleavage_counts_t noisy;
  cleavage_counts_clear(&noisy);
  add_n(&noisy, HPAI_NT, 100);
  char nt[CLEAVAGE_SITE_NT + 1];
  memcpy(nt, "TTTATAGAGGGAGGATGGCAGGGA", sizeof(nt));
  for (size_t i = 0; i < 100; i++) {
    for (size_t j = 0; j < 4; j++)
      nt[j] = "ACGT"[(i >> (2 * j)) & 3];
    for (size_t w = 0; w < CLEAVAGE_N_CONTROLS; w++)
      cleavage_counts_add_control(&noisy, w, nt);
  }
  cleavage_result_of(&noisy, CLEAVAGE_CHEMISTRY_R10, &r);
  test_assert(r.control_support < CLEAVAGE_MIN_CONTROL_SUPPORT);
  test_assert(r.control_weak);
  test_eq_size(r.verdict, CLEAVAGE_CONFIRM);
  /* the site itself is still reported, with its own support */
  test_eq_str(r.aa, "QRERRRKR");
  test_assert(r.support > 0.99);
}

/* the fusion peptide: the anchor as the reads carry it */
static void test_fusion(void) {
  cleavage_counts_t c;
  cleavage_counts_clear(&c);
  cleavage_counts_set_reference(&c, "LREKRRKR", "GLFGAIAG");
  add_n(&c, HPAI_NT, 100);
  cleavage_result_t r;
  /* no read read at the anchor yet */
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_str(r.reference_fusion, "GLFGAIAG");
  test_eq_str(r.fusion_aa, "");
  test_eq_size(r.fusion_reads, 0);
  /* 80 reads carry the peptide of the reference, 20 another */
  for (size_t i = 0; i < 80; i++)
    cleavage_counts_add_fusion(&c, "GGACTGTTTGGGGCTATAGCAGGC");
  for (size_t i = 0; i < 20; i++)
    cleavage_counts_add_fusion(&c, "GGAATTTTTGGGGCTATAGCAGGC");
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_str(r.fusion_aa, "GLFGAIAG");
  test_eq_size(r.fusion_reads, 100);
  test_assert(r.fusion_support > 0.79 && r.fusion_support < 0.81);
  /* the other peptide is the one of the H3 of horses and pigs */
  cleavage_counts_t h3;
  cleavage_counts_clear(&h3);
  for (size_t i = 0; i < 50; i++)
    cleavage_counts_add_fusion(&h3, "GGAATTTTTGGGGCTATAGCAGGC");
  add_n(&h3, HPAI_NT, 50);
  cleavage_result_of(&h3, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_str(r.fusion_aa, "GIFGAIAG");
  test_assert(r.fusion_support > 0.99);
}

static void test_readings(void) {
  cleavage_counts_t c;
  cleavage_result_t r;
  /* a low pathogenic site: one basic residue before HA2 */
  cleavage_counts_clear(&c);
  add_n(&c, LPAI_NT, 50);
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_str(r.aa, "NVPARSSR");
  test_eq_size(r.reading, CLEAVAGE_MONOBASIC);
  test_eq_size(r.verdict, CLEAVAGE_FIRM);
  /* a motif with a stop codon cannot be read */
  cleavage_counts_clear(&c);
  add_n(&c, "CAAAGAGAGTAGAGAAGAAAAAGA", 50);
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_size(r.reading, CLEAVAGE_UNREADABLE);
  /* and one with an ambiguous base */
  cleavage_counts_clear(&c);
  add_n(&c, "CAAAGAGAGAGAAGAAGAAAAANA", 50);
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_size(r.reading, CLEAVAGE_UNREADABLE);
}

static void test_overflow(void) {
  /* more distinct motifs than the table holds: the reads are counted,
     the support is a lower bound, and the frequent motif is kept */
  cleavage_counts_t c;
  cleavage_counts_clear(&c);
  add_n(&c, HPAI_NT, 500);
  /* the other motifs are built on a different template, so that none of
     them is the frequent one */
  char nt[CLEAVAGE_SITE_NT + 1];
  memcpy(nt, LPAI_NT, sizeof(nt));
  for (size_t i = 0; i < CLEAVAGE_MAX_MOTIFS + 100; i++) {
    for (size_t j = 0; j < 6; j++)
      nt[j] = "ACGT"[(i >> (2 * j)) & 3];
    add_n(&c, nt, 1);
  }
  test_eq_size(c.n_motifs, CLEAVAGE_MAX_MOTIFS);
  test_assert(c.overflow > 0);
  cleavage_result_t r;
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_str(r.nt, HPAI_NT);
  test_eq_size(r.motif_reads, 500);
}

static void test_save_load(void) {
  cleavage_counts_t c;
  cleavage_counts_clear(&c);
  add_n(&c, HPAI_NT, 30);
  add_n(&c, LPAI_NT, 7);
  c.overflow = 3;
  char *buf = NULL;
  size_t size = 0;
  FILE *f = open_memstream(&buf, &size);
  test_assert(f != NULL);
  cleavage_counts_save(&c, f);
  fclose(f);
  cleavage_counts_t loaded;
  cleavage_counts_clear(&loaded);
  for (char *line = strtok(buf, "\n"); line != NULL;
       line = strtok(NULL, "\n"))
    test_assert(cleavage_counts_load_line(&loaded, line));
  test_eq_size(loaded.searched, 37);
  test_eq_size(loaded.overflow, 3);
  test_eq_size(cleavage_counts_reads(&loaded), 40);
  test_eq_size(loaded.n_motifs, 2);
  free(buf);
  /* a line of another kind is not one of these */
  char other[] = "reads 1 2 3";
  test_assert(!cleavage_counts_load_line(&loaded, other));
}

static void test_chemistry(void) {
  test_eq_size(cleavage_chemistry(true, NULL), CLEAVAGE_CHEMISTRY_ILLUMINA);
  test_eq_size(cleavage_chemistry(true, "dna_r9.4.1_e8_hac@v3.3"),
               CLEAVAGE_CHEMISTRY_ILLUMINA);
  test_eq_size(cleavage_chemistry(false, NULL), CLEAVAGE_CHEMISTRY_UNKNOWN);
  test_eq_size(cleavage_chemistry(false, ""), CLEAVAGE_CHEMISTRY_UNKNOWN);
  test_eq_size(cleavage_chemistry(false, "dna_r10.4.1_e8.2_400bps_hac@v4.2.0"),
               CLEAVAGE_CHEMISTRY_R10);
  test_eq_size(cleavage_chemistry(false, "dna_r9.4.1_e8_hac@v3.3"),
               CLEAVAGE_CHEMISTRY_R9);
  test_eq_size(cleavage_chemistry(false, "something_else"),
               CLEAVAGE_CHEMISTRY_UNKNOWN);
  test_eq_str(cleavage_chemistry_name(CLEAVAGE_CHEMISTRY_R10), "R10");
  test_eq_str(cleavage_reading_name(CLEAVAGE_MULTIBASIC), "multibasic");
  test_eq_str(cleavage_verdict_name(CLEAVAGE_CONFIRM), "confirm");
}

void cleavage_counts_test(void) {
  test_counts();
  test_supported_readings();
  test_control();
  test_fusion();
  test_readings();
  test_overflow();
  test_save_load();
  test_chemistry();
}
