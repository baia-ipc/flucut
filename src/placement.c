/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <string.h>
#include "fasta.h"
#include "placement.h"
#include "primers.h"
#include "test.h"

/* the span of a read of <len> bases whose probe, taken at <offset>,
   was located as <hit> (see reference_hit_t: for a reverse hit,
   <ref_start> is the position on the plus strand of the first base of
   the probe, and the read runs towards the start of the reference) */
static void span_of(const reference_hit_t *hit, long offset, long len,
    placement_t *out) {
  out->found = true;
  out->ref_index = hit->ref_index;
  out->ref_len = hit->ref_len;
  out->reverse = hit->reverse;
  if (hit->reverse) {
    out->to = hit->ref_start + offset + 1;
    out->from = out->to - len;
  } else {
    out->from = hit->ref_start - offset;
    out->to = out->from + len;
  }
}

void placement_locate(reference_set_t *refs, const char *seq, size_t len,
    placement_t *out) {
  assert(refs != NULL);
  assert(seq != NULL || len == 0);
  assert(out != NULL);
  memset(out, 0, sizeof(placement_t));
  if (len < REFERENCE_MIN_PROBE_LEN)
    return;
  size_t probe_len = reference_set_probe_len(refs);
  if (probe_len > len)
    probe_len = len;
  /* the probes: in the middle of the read, then at a quarter and at
     three quarters of it (the same probe, in a short read, is tried
     once) */
  size_t room = len - probe_len;
  const size_t offsets[PLACEMENT_N_PROBES] = {room / 2, room / 4,
                                              room - room / 4};
  for (size_t i = 0; i < PLACEMENT_N_PROBES; i++) {
    if (i > 0 && offsets[i] == offsets[i - 1])
      continue;
    reference_hit_t hit;
    reference_set_locate_seeded(refs, seq + offsets[i], probe_len, &hit);
    if (hit.mapped) {
      span_of(&hit, (long)offsets[i], (long)len, out);
      return;
    }
  }
}

/* --- unit tests --- */

void placement_test(void) {
  enum { REF_LEN = 1000, PROBE = 32, READ = 300 };
  unsigned long state = 7;
  char ref0[REF_LEN + 1], ref1[REF_LEN + 1], read[READ + 1];
  random_seq(ref0, REF_LEN, &state);
  random_seq(ref1, REF_LEN, &state);
  fasta_record_t recs[2] = {{"r0", "", ref0, REF_LEN},
                            {"r1", "", ref1, REF_LEN}};
  fasta_t fa = {recs, 2};
  reference_index_t *ix = reference_index_new(&fa);
  reference_set_t *rs = reference_set_new(ix, PROBE, 0.2);
  placement_t p;

  /* a piece of the interior of a reference */
  memcpy(read, ref1 + 400, READ);
  placement_locate(rs, read, READ, &p);
  test_assert(p.found);
  test_eq_size(p.ref_index, 1);
  test_assert(!p.reverse);
  test_eq_long(p.from, 400);
  test_eq_long(p.to, 700);
  test_eq_long(p.ref_len, REF_LEN);

  /* the same piece, read from the other strand: the same span */
  primers_reverse_complement(ref1 + 400, READ, read);
  placement_locate(rs, read, READ, &p);
  test_assert(p.found);
  test_eq_size(p.ref_index, 1);
  test_assert(p.reverse);
  test_eq_long(p.from, 400);
  test_eq_long(p.to, 700);

  /* the middle of the read is not in the references: the read is
     placed by the probe at a quarter of it */
  memcpy(read, ref0 + 100, READ);
  random_seq(read + 120, 60, &state);
  read[180] = ref0[280];   /* (random_seq() ends the sequence) */
  placement_locate(rs, read, READ, &p);
  test_assert(p.found);
  test_eq_size(p.ref_index, 0);
  test_eq_long(p.from, 100);
  test_eq_long(p.to, 400);

  /* a read longer than what is left of the reference: the span is not
     limited to the reference */
  memcpy(read, ref0 + 900, 100);
  random_seq(read + 100, READ - 100, &state);
  placement_locate(rs, read, 110, &p);
  test_assert(p.found);
  test_eq_long(p.from, 900);
  test_eq_long(p.to, 1010);

  /* a read shorter than the probe is placed whole; too short, not */
  placement_locate(rs, ref0 + 500, 20, &p);
  test_assert(p.found);
  test_eq_long(p.from, 500);
  test_eq_long(p.to, 520);
  placement_locate(rs, ref0 + 500, REFERENCE_MIN_PROBE_LEN - 1, &p);
  test_assert(!p.found);

  /* a sequence which is not in the references */
  random_seq(read, READ, &state);
  placement_locate(rs, read, READ, &p);
  test_assert(!p.found);

  reference_set_delete(rs);
  reference_index_delete(ix);
}
