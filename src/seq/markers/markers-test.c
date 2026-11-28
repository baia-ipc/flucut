/* (c) 2026 Giorgio Gonnella */

#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "marker-counts.h"
#include "marker-read.h"
#include "markers-test.h"
#include "markers.h"
#include "reference-index.h"
#include "sequtil.h"
#include "test.h"

/* the synthetic references of the tests: the anchor of a marker in a
   little sequence of its segment, with what the marker reads beside it.
   The codons are the plainest ones of each residue */

/* PB2: MDVNPTLLFLKV, then 627 (E), then the anchor QSRMQFSS, LTV and a
   tail */
#define PB2_HEAD "ATGGATGTTAATCCTACTCTTCTTTTTCTTAAAGTT"
#define PB2_ANCHOR "CAATCTCGTATGCAATTTTCTTCT"
#define PB2_TAIL "CTTACTGTTAATGCTCCTAATTATAAAAAACTTTCTATGAATTCT"
#define PB2_627E PB2_HEAD "GAA" PB2_ANCHOR PB2_TAIL
#define PB2_627K PB2_HEAD "AAA" PB2_ANCHOR PB2_TAIL

/* NS1: the anchor DGRPSLPP and the C-terminus EKK* after it */
#define NS1_ANCHOR "GATGGTCGTCCTTCTCTTCCTCCT"
#define NS1_SEQ PB2_HEAD NS1_ANCHOR "GAAAAAAAATAA" \
                "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT"

/* NA: the anchor at its start, a stalk of ten residues and the second
   anchor HSNGTVKDR in the head */
#define NA_START "ATGAATCCTAATCAAAAAATTACT"
#define NA_STALK "GTTTCTCTTACTATTGCTACTATTTGTTTT"
#define NA_HEAD "CATTCTAATGGTACTGTTAAAGATCGT"
#define NA_SEQ NA_START NA_STALK NA_HEAD "GCTGCTGCTGCTGCT"

static const char *const SEQS[] = { PB2_627E, NS1_SEQ, NA_SEQ };
static const char *const DESCS[] = {
  "Influenza A virus (A/test/2020(H5N1)) segment 1 polymerase PB2 (PB2) "
  "gene, complete cds",
  "Influenza A virus (A/test/2020(H5N1)) segment 8 nonstructural protein "
  "(NS) gene, complete cds",
  "Influenza A virus (A/test/2020(H5N1)) segment 6 neuraminidase (NA) "
  "gene, complete cds",
};
#define N_SEQS (sizeof(SEQS) / sizeof(SEQS[0]))

static fasta_t* test_fasta(void) {
  fasta_t *fa;
  XMALLOC(fa, sizeof(fasta_t));
  XMALLOC(fa->records, sizeof(fasta_record_t) * N_SEQS);
  fa->n_records = N_SEQS;
  for (size_t i = 0; i < N_SEQS; i++) {
    fa->records[i].name = (char*)"ref";
    fa->records[i].description = (char*)DESCS[i];
    fa->records[i].seq = (char*)SEQS[i];
    fa->records[i].seq_len = strlen(SEQS[i]);
  }
  return fa;
}

/* --- the values read from a window --- */

static void test_values(void) {
  char v[MARKER_VALUE_MAX + 1];
  /* the residues of a window */
  test_assert(marker_value_of(MARKER_RESIDUES, "GAA", 3, v));
  test_eq_str(v, "E");
  test_assert(marker_value_of(MARKER_RESIDUES, "CAATCTGGT", 9, v));
  test_eq_str(v, "QSG");
  /* a codon with an ambiguous base is an X, and is a value like any
     other: the reads carrying it are not the ones agreeing */
  test_assert(marker_value_of(MARKER_RESIDUES, "GNA", 3, v));
  test_eq_str(v, "X");
  /* the residues up to the stop codon, the stop included */
  test_assert(marker_value_of(MARKER_TO_STOP, "GAAAAAAAATAAGCTGCT", 18, v));
  test_eq_str(v, "EKK*");
  /* a window without a stop codon does not say where the protein ends */
  test_assert(!marker_value_of(MARKER_TO_STOP, "GAAAAAAAAGCTGCTGCT", 18, v));
  /* a protein which ends farther into the window gives its last
     residues, so that the value stays a few codons wide */
  test_assert(marker_value_of(MARKER_TO_STOP,
      "AAACAAAAACGTAAAATGGCTCGTACTATTGAATCTGAAGTTTAAGCT", 48, v));
  test_eq_str(v, "ESEV*");
  /* the codons before the stop codon */
  test_assert(marker_value_of(MARKER_FRAME_LENGTH, "GAAAAAAAATAAGCT", 15, v));
  test_eq_str(v, "3");
  test_assert(!marker_value_of(MARKER_FRAME_LENGTH, "GAAAAAAAAGCTGCT", 15, v));
  marker_span_value(327, v);
  test_eq_str(v, "327");

  /* a value deserves attention when it has one of the residues of the
     table at the position where the table has it */
  const marker_def_t *pb2 = marker_def(MARKER_PB2_627);
  test_assert(marker_value_attention(pb2, "K", ""));
  test_assert(!marker_value_attention(pb2, "E", ""));
  /* the haemagglutinin is read at three positions, of which the first
     and the third are the ones watched */
  const marker_def_t *ha = marker_def(MARKER_HA_226_228);
  test_assert(marker_value_attention(ha, "LSS", ""));
  test_assert(marker_value_attention(ha, "QSS", ""));
  test_assert(marker_value_attention(ha, "LSG", ""));
  test_assert(!marker_value_attention(ha, "QSG", ""));
}

/* --- the anchors in the references --- */

static void test_anchors(void) {
  fasta_t *fa = test_fasta();
  reference_index_t *ix = reference_index_new(fa);
  marker_anchors_t *anchors = marker_anchors_new(fa, ix, 0.2, NULL);
  test_assert(anchors != NULL);
  /* the three markers whose anchors are in these three references */
  test_eq_size(marker_anchors_n(anchors), 3);
  const marker_reference_t *chosen[MARKER_N_MARKERS];
  marker_anchors_chosen(anchors, NULL, chosen);
  test_assert(chosen[MARKER_PB2_627] != NULL);
  test_assert(chosen[MARKER_NS1_TAIL] != NULL);
  test_assert(chosen[MARKER_NA_STALK] != NULL);
  /* the segments without a reference have none */
  test_assert(chosen[MARKER_HA_226_228] == NULL);
  test_assert(chosen[MARKER_PA_38] == NULL);
  /* what each reference itself carries */
  test_eq_str(chosen[MARKER_PB2_627]->value, "E");
  test_eq_str(chosen[MARKER_PB2_627]->anchor_aa, "QSRMQFSS");
  test_eq_str(chosen[MARKER_NS1_TAIL]->value, "EKK*");
  /* the stalk: the bases from the start of NA to its second anchor */
  test_eq_str(chosen[MARKER_NA_STALK]->value, "54");
  /* only the references allowed are used */
  bool allowed[N_SEQS] = { false, true, true };
  marker_anchors_chosen(anchors, allowed, chosen);
  test_assert(chosen[MARKER_PB2_627] == NULL);
  test_assert(chosen[MARKER_NS1_TAIL] != NULL);
  marker_anchors_delete(anchors);
  reference_index_delete(ix);
  free(fa->records);
  free(fa);
}

/* --- the markers read from the reads --- */

/* count <n> reads of <seq>, and of its reverse complement, located on
   the reference <ref>, with the junction <junction> (-1: none) */
static void read_n_at(const marker_reference_t *const *chosen, long ref,
    long junction, const char *seq, size_t n, marker_counts_t *counts) {
  char rc[512];
  size_t len = strlen(seq);
  marker_breaks_t breaks = {&junction, junction < 0 ? 0 : 1, NULL, 0};
  reverse_complement(seq, len, rc);
  for (size_t i = 0; i < n; i++) {
    marker_read_sites(chosen, ref, &breaks, seq, len, 0.2, counts);
    marker_read_sites(chosen, ref, &breaks, rc, len, 0.2, counts);
  }
}

static void read_n(const marker_reference_t *const *chosen, long ref,
    const char *seq, size_t n, marker_counts_t *counts) {
  read_n_at(chosen, ref, -1, seq, n, counts);
}

static void test_reads(void) {
  fasta_t *fa = test_fasta();
  reference_index_t *ix = reference_index_new(fa);
  marker_anchors_t *anchors = marker_anchors_new(fa, ix, 0.2, NULL);
  const marker_reference_t *chosen[MARKER_N_MARKERS];
  marker_anchors_chosen(anchors, NULL, chosen);
  marker_counts_t c;
  marker_counts_clear(&c);
  /* 45 reads carrying the residue of the reference and 5 carrying the
     other one, each read in both orientations */
  read_n(chosen, 0, PB2_627E, 45, &c);
  read_n(chosen, 0, PB2_627K, 5, &c);
  test_eq_size(marker_counts_reads(&c, MARKER_PB2_627), 100);
  marker_result_t r;
  marker_result_of(&c, MARKER_PB2_627, &r);
  test_assert(r.found);
  test_eq_size(r.searched, 100);
  test_eq_size(r.reads, 100);
  test_eq_str(r.value, "E");
  test_eq_size(r.value_reads, 90);
  test_eq_size(r.distinct, 2);
  test_eq_str(r.reference_value, "E");
  test_assert(!r.differs_from_reference);
  test_eq_size(r.verdict, MARKER_FIRM);
  /* the reads agree on the anchor, so the minority is not noise */
  test_eq_size(r.anchor_reads, 100);
  test_assert(!r.anchor_weak);
  test_assert(r.anchor_support > 0.99);
  /* the residue a tenth of the reads carry is reported, with its share,
     which is what reading the position off a consensus would lose */
  test_eq_size(r.supported, 2);
  test_eq_size(r.n_top, 2);
  test_eq_str(r.top[0].value, "E");
  test_eq_str(r.top[1].value, "K");
  test_eq_size(r.top[1].reads, 10);
  test_assert(r.top[1].share > 0.09 && r.top[1].share < 0.11);
  test_assert(r.top[1].attention);
  test_assert(!r.top[0].attention);
  test_assert(r.attention);
  test_assert(r.attention_share > 0.09 && r.attention_share < 0.11);

  /* a read located on another reference is not searched for a marker
     of this one */
  marker_counts_clear(&c);
  read_n(chosen, 2, PB2_627E, 10, &c);
  marker_result_of(&c, MARKER_PB2_627, &r);
  test_eq_size(r.searched, 0);
  /* a read which is not located is searched for every marker */
  marker_counts_clear(&c);
  read_n(chosen, -1, PB2_627E, 5, &c);
  marker_result_of(&c, MARKER_PB2_627, &r);
  test_eq_size(r.searched, 10);
  test_eq_size(r.reads, 10);
  /* too few reads for a firm value */
  test_eq_size(r.verdict, MARKER_TOO_FEW_READS);

  /* reads which agree on no value: none is reported as the sample's */
  marker_counts_clear(&c);
  for (size_t i = 0; i < 100; i++) {
    char value[4] = {"ACDEFGHIKLMNPQRSTVWY"[i % 20], '\0', '\0', '\0'};
    marker_counts_searched(&c, MARKER_PB2_627);
    marker_counts_add(&c, MARKER_PB2_627, value);
    marker_counts_add_anchor(&c, MARKER_PB2_627, "QSRMQFSS");
  }
  marker_result_of(&c, MARKER_PB2_627, &r);
  test_eq_size(r.verdict, MARKER_UNRESOLVED);

  /* a read with an internal deletion gives the marker only where its
     anchor and its window lie on the same side of the junction: the
     bases before the junction come from another part of the segment */
  marker_counts_clear(&c);
  size_t pb2_len = strlen(PB2_627E);
  long at = (long)strlen(PB2_HEAD) + 3;   /* the first base of the anchor */
  marker_breaks_t at_anchor = {&at, 1, NULL, 0};
  for (size_t i = 0; i < 30; i++)
    marker_read_sites(chosen, 0, &at_anchor, PB2_627E, pb2_len, 0.2, &c);
  marker_result_of(&c, MARKER_PB2_627, &r);
  test_eq_size(r.searched, 30);
  /* the window lies before the junction and the anchor after it: no
     value, though the anchor was found */
  test_eq_size(r.reads, 0);
  test_assert(r.anchor_reads > 0);
  /* a junction past the whole site: the value is read */
  marker_counts_clear(&c);
  long past = (long)pb2_len;
  marker_breaks_t beyond = {&past, 1, NULL, 0};
  for (size_t i = 0; i < 30; i++)
    marker_read_sites(chosen, 0, &beyond, PB2_627E, pb2_len, 0.2, &c);
  marker_result_of(&c, MARKER_PB2_627, &r);
  test_eq_size(r.reads, 30);
  test_eq_str(r.value, "E");

  /* a primer which matched inside the read: its bases are the primer
     and not the template, so a window which overlaps it gives nothing */
  marker_counts_clear(&c);
  long at_primer = (long)strlen(PB2_HEAD);
  marker_blocked_t primer = {at_primer, at_primer + 2};
  marker_breaks_t inside = {NULL, 0, &primer, 1};
  for (size_t i = 0; i < 30; i++)
    marker_read_sites(chosen, 0, &inside, PB2_627E, pb2_len, 0.2, &c);
  marker_result_of(&c, MARKER_PB2_627, &r);
  test_eq_size(r.searched, 30);
  test_eq_size(r.reads, 0);
  test_assert(r.anchor_reads > 0);

  /* the C-terminus of a protein, and the length of a stretch, read in
     both orientations */
  marker_counts_clear(&c);
  read_n(chosen, 1, NS1_SEQ, 20, &c);
  marker_result_of(&c, MARKER_NS1_TAIL, &r);
  test_eq_str(r.value, "EKK*");
  test_eq_size(r.reads, 40);
  marker_counts_clear(&c);
  read_n(chosen, 2, NA_SEQ, 20, &c);
  marker_result_of(&c, MARKER_NA_STALK, &r);
  test_eq_str(r.value, "54");
  test_eq_size(r.reads, 40);
  test_assert(!r.differs_from_reference);

  marker_anchors_delete(anchors);
  reference_index_delete(ix);
  free(fa->records);
  free(fa);
}

/* --- merged, saved and loaded --- */

static void test_save_load(void) {
  marker_counts_t a, b;
  marker_counts_clear(&a);
  marker_counts_clear(&b);
  marker_counts_set_reference(&a, MARKER_PB2_627, "E", "QSRMQFSS");
  for (size_t i = 0; i < 30; i++) {
    marker_counts_searched(&a, MARKER_PB2_627);
    marker_counts_add(&a, MARKER_PB2_627, i < 25 ? "E" : "K");
    marker_counts_add_anchor(&a, MARKER_PB2_627, "QSRMQFSS");
  }
  for (size_t i = 0; i < 10; i++) {
    marker_counts_searched(&b, MARKER_PB2_627);
    marker_counts_add(&b, MARKER_PB2_627, "E");
    marker_counts_add_anchor(&b, MARKER_PB2_627, "QSRMQFSS");
  }
  marker_counts_merge(&a, &b);
  marker_result_t r;
  marker_result_of(&a, MARKER_PB2_627, &r);
  test_eq_size(r.reads, 40);
  test_eq_size(r.value_reads, 35);

  /* saved and read back: the same counts */
  char *buf = NULL;
  size_t size = 0;
  FILE *file = open_memstream(&buf, &size);
  marker_counts_save(&a, file);
  fclose(file);
  marker_counts_t loaded;
  marker_counts_clear(&loaded);
  for (char *line = strtok(buf, "\n"); line != NULL;
       line = strtok(NULL, "\n"))
    test_assert(marker_counts_load_line(&loaded, line));
  marker_result_t r2;
  marker_result_of(&loaded, MARKER_PB2_627, &r2);
  test_eq_size(r2.searched, r.searched);
  test_eq_size(r2.reads, r.reads);
  test_eq_size(r2.value_reads, r.value_reads);
  test_eq_str(r2.value, r.value);
  test_eq_str(r2.reference_value, "E");
  test_assert(r2.anchor_support > 0.99);
  /* a line of another kind is not one of these */
  char other[] = "cleavage-reads 10 0";
  test_assert(!marker_counts_load_line(&loaded, other));
  free(buf);
}

/* the reads whose value finds no slot are counted as reads, which can
   only lower the share reported */
static void test_overflow(void) {
  marker_counts_t c;
  marker_counts_clear(&c);
  char value[MARKER_VALUE_MAX + 1];
  for (size_t i = 0; i < MARKER_VALUES + 50; i++) {
    marker_counts_searched(&c, MARKER_PB1_F2);
    snprintf(value, sizeof(value), "%zu", i);
    marker_counts_add(&c, MARKER_PB1_F2, value);
  }
  for (size_t i = 0; i < 100; i++)
    marker_counts_add(&c, MARKER_PB1_F2, "90");
  marker_result_t r;
  marker_result_of(&c, MARKER_PB1_F2, &r);
  test_eq_size(r.distinct, MARKER_VALUES);
  test_eq_size(r.reads, MARKER_VALUES + 150);
  test_eq_str(r.value, "90");
}

void markers_test(void) {
  test_values();
  test_anchors();
  test_reads();
  test_save_load();
  test_overflow();
}
