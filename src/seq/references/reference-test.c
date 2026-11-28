/* (c) 2026 Giorgio Gonnella */

/*
   Unit tests of the references: the group of a reference
   (reference-subtype.h), and the location of a probe
   (reference.h), with the k-mer index (reference-index.h) or without.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "fasta.h"
#include "primers.h"
#include "reference.h"
#include "reference-index.h"
#include "reference-subtype.h"
#include "reference-test.h"
#include "test.h"

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

/* the tests of the location of a probe, with the index or exhaustive
   (the same results; <ref0>, <ref1>: the references of <ix>, see
   reference_test()) */
static void test_locate(reference_index_t *ix, bool exhaustive,
    const char *ref0, const char *ref1, unsigned long *state) {
  enum { REF_LEN = 1000, PROBE = 30 };
  char probe[64];
  reference_set_t *rs = reference_set_new(ix, PROBE, 0.2);
  reference_set_set_exhaustive(rs, exhaustive);
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
  test_random_dna_seq(probe, PROBE, state);
  reference_set_locate(rs, probe, PROBE, &hit);
  test_assert(!hit.mapped);
  test_eq_long(hit.ref_start, -1);

  /* a probe too short to be located */
  take_probe(ref1, 500, REFERENCE_MIN_PROBE_LEN - 1, false, probe);
  reference_set_locate(rs, probe, REFERENCE_MIN_PROBE_LEN - 1, &hit);
  test_assert(!hit.mapped);

  reference_set_delete(rs);
}

void reference_test(void) {
  /* the group of a reference: the subtype for influenza A, the lineage
     for influenza B, and "B" where no lineage is named */
  char group[32];
  test_assert(reference_description_subtype("Influenza A virus "
      "(A/Cambodia/1/2023(H5N1)) segment 4", group, sizeof(group)));
  test_eq_str(group, "H5N1");
  test_assert(reference_description_subtype("Influenza B virus "
      "(B/Brisbane/60/2008) segment 4", group, sizeof(group)));
  test_eq_str(group, "B");
  test_assert(reference_description_subtype("Influenza B virus "
      "(B/Phuket/3073/2013) Yamagata lineage, segment 4", group,
      sizeof(group)));
  test_eq_str(group, "Yamagata");
  test_assert(reference_description_subtype("Influenza B virus, "
      "lineage Victoria, segment 6", group, sizeof(group)));
  test_eq_str(group, "Victoria");
  /* the place of a strain name is not a lineage */
  test_assert(reference_description_subtype("Influenza B virus "
      "(B/Victoria/504/2000) segment 4", group, sizeof(group)));
  test_eq_str(group, "B");
  /* a record of influenza A without a subtype has no group */
  test_assert(!reference_description_subtype("Influenza A virus "
      "(A/Hong Kong/1/68) segment 4", group, sizeof(group)));
  enum { REF_LEN = 1000, PROBE = 30 };
  unsigned long state = 42;
  char ref0[REF_LEN + 1], ref1[REF_LEN + 1], probe[64];
  test_random_dna_seq(ref0, REF_LEN, &state);
  test_random_dna_seq(ref1, REF_LEN, &state);
  /* a repeat in ref1: the first 40 bases are also at 600 */
  memcpy(ref1 + 600, ref1, 40);
  /* a degenerate base in ref0 */
  ref0[300] = 'R';
  char *seqs[] = {ref0, ref1};
  fasta_t *fa = test_fasta(seqs, 2);
  reference_index_t *ix = reference_index_new(fa);
  test_assert(reference_index_fasta(ix) == fa);
  test_locate(ix, false, ref0, ref1, &state);
  test_locate(ix, true, ref0, ref1, &state);

  /* a probe with an error in each k-mer: the k-mers find nothing, the
     exhaustive search which follows finds it */
  take_probe(ref1, 500, PROBE, false, probe);
  for (size_t i = REFERENCE_KMER_LEN - 1; i < PROBE;
       i += REFERENCE_KMER_LEN - 1)
    probe[i] = probe[i] == 'A' ? 'C' : 'A';
  reference_set_t *rs = reference_set_new(ix, PROBE, 0.2);
  reference_hit_t hit;
  reference_set_locate(rs, probe, PROBE, &hit);
  test_assert(hit.mapped);
  test_eq_long(hit.ref_start, 500);
  test_eq_long(hit.n_errors, 3);
  /* by the k-mers alone it is not found, whatever the settings */
  reference_set_locate_seeded(rs, probe, PROBE, &hit);
  test_assert(!hit.mapped);
  reference_set_set_exhaustive(rs, true);
  reference_set_locate_seeded(rs, probe, PROBE, &hit);
  test_assert(!hit.mapped);
  reference_set_set_exhaustive(rs, false);
  /* with the seeds only, it is not found */
  reference_set_set_seeds_only(rs, true);
  reference_set_locate(rs, probe, PROBE, &hit);
  test_assert(!hit.mapped);
  /* ... while a probe with an exact k-mer is */
  take_probe(ref1, 500, PROBE, false, probe);
  probe[3] = probe[3] == 'A' ? 'C' : 'A';
  reference_set_locate(rs, probe, PROBE, &hit);
  test_assert(hit.mapped);
  test_eq_long(hit.ref_start, 500);
  reference_set_set_seeds_only(rs, false);
  reference_set_delete(rs);
  reference_index_delete(ix);

  /* the references of one segment only: a probe of another segment is
     not located, by the k-mers or exhaustively */
  char ha[] = "segment 4 (HA)", na[] = "segment 6 (NA)";
  fa->records[0].description = ha;
  fa->records[1].description = na;
  ix = reference_index_new(fa);
  rs = reference_set_new(ix, PROBE, 0.2);
  take_probe(ref1, 500, PROBE, false, probe);
  reference_set_restrict_segment(rs, 4);
  reference_set_locate(rs, probe, PROBE, &hit);
  test_assert(!hit.mapped);
  reference_set_set_exhaustive(rs, true);
  reference_set_locate(rs, probe, PROBE, &hit);
  test_assert(!hit.mapped);
  reference_set_set_exhaustive(rs, false);
  reference_set_restrict_segment(rs, 6);
  reference_set_locate(rs, probe, PROBE, &hit);
  test_assert(hit.mapped);
  test_eq_size(hit.ref_index, 1);
  test_eq_long(hit.ref_start, 500);
  reference_set_restrict_segment(rs, 0);
  reference_set_locate(rs, probe, PROBE, &hit);
  test_assert(hit.mapped);
  reference_set_delete(rs);
  reference_index_delete(ix);
  free(fa->records);
  free(fa);
}
