/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "chimera.h"
#include "fasta.h"
#include "junction.h"
#include "primers.h"
#include "terminus-matches.h"
#include "test.h"

struct chimera_finder_t {
  terminus_matches_t *termini;
  long window;
  reference_set_t *refs;
  /* the bases next to a match, on the plus strand of the reference */
  char *bases;
  size_t bases_cap;
};

chimera_finder_t* chimera_finder_new(const scheme_t *scheme,
    double max_errors, long window) {
  assert(scheme != NULL);
  assert(window >= 0);
  chimera_finder_t *self;
  XMALLOC(self, sizeof(chimera_finder_t));
  self->termini = terminus_matches_new(scheme, max_errors);
  self->window = window;
  self->refs = NULL;
  self->bases = NULL;
  self->bases_cap = 0;
  return self;
}

void chimera_finder_restrict_virus(chimera_finder_t *self,
    virus_set_t types) {
  assert(self != NULL);
  terminus_matches_restrict_virus(self->termini, types);
}

void chimera_finder_use_references(chimera_finder_t *self,
    reference_set_t *refs) {
  assert(self != NULL);
  assert(refs != NULL);
  self->refs = refs;
}

void chimera_finder_delete(chimera_finder_t *self) {
  if (self == NULL)
    return;
  terminus_matches_delete(self->termini);
  free(self->bases);
  free(self);
}

/* locate in the references the bases of <seq> (length <len>) from
   <from> to <to> - 1, limited to the read (at most a probe of them) */
static void locate(chimera_finder_t *self, const char *seq, size_t len,
    long from, long to, reference_hit_t *hit) {
  hit->mapped = false;
  if (from < 0)
    from = 0;
  if (to > (long)len)
    to = (long)len;
  if (to - from >= REFERENCE_MIN_PROBE_LEN)
    reference_set_locate(self->refs, seq + from, (size_t)(to - from), hit);
}

/*
   Whether two probes of a read, <a> and then <b>, whose first bases
   are <in_read> bases apart, are located on the same reference and
   strand, at the distance they have in the read, give or take
   <tolerance>, or farther apart in the same order: the read then lacks
   bases of the reference between them (an internal deletion, see
   amplicon.h), and is one template all the same.
*/
static bool is_continuous(const reference_hit_t *a,
    const reference_hit_t *b, long in_read, long tolerance) {
  if (!a->mapped || !b->mapped || a->ref_index != b->ref_index ||
      a->reverse != b->reverse)
    return false;
  long in_ref = b->reverse ? a->ref_start - b->ref_start :
                b->ref_start - a->ref_start;
  return in_ref - in_read >= -tolerance;
}

/*
   Whether the bases on the other side of the match <m> continue the
   reference of <inner> (the probe of the amplicon side, which mapped):
   they are aligned to that reference, from where the probe lies,
   outwards across the match (CHIMERA_OUTER_PROBES probes of them, the
   bases of the match included, which belong to the segment if the read
   runs through them). The read then runs through the match.

   The alignment tolerates the differences of a sample from its
   reference and the errors of a read (CHIMERA_MAX_DIFFERENCE), which a
   search of those bases in the references does not: a probe of a sample
   far from every reference is often not located at all, and the absence
   of a location is no evidence that the read stops there. Bases of
   another amplicon, which do not belong next to the probe, differ far
   more.
*/
static bool continues_reference(chimera_finder_t *self, const char *seq,
    size_t len, const terminus_match_t *m, const reference_hit_t *inner) {
  long probe = (long)reference_set_probe_len(self->refs);
  long match_len = m->to - m->from + 1;
  long from, to, anchor;
  bool backwards;
  if (m->start) {
    /* the amplicon side is after the match: the bases before it, with
       those of the match, end where the probe begins */
    to = m->to + 1;
    from = to - CHIMERA_OUTER_PROBES * probe - match_len;
    if (from < 0)
      from = 0;
    anchor = inner->reverse ? inner->ref_start + 1 : inner->ref_start - 1;
    backwards = !inner->reverse;
  } else {
    /* the amplicon side is before the match: the bases from the match
       on begin where the probe ends */
    from = m->from;
    to = from + CHIMERA_OUTER_PROBES * probe + match_len;
    if (to > (long)len)
      to = (long)len;
    anchor = inner->reverse ? inner->ref_start - probe :
             inner->ref_start + probe;
    backwards = inner->reverse;
  }
  long n = to - from;
  if (n < REFERENCE_MIN_PROBE_LEN)
    return false;   /* too few bases to tell */
  const char *plus = seq + from;
  if (inner->reverse) {
    XCHECKSPACEINC(self->bases, &self->bases_cap, (size_t)n + 1, 256);
    primers_reverse_complement(seq + from, (size_t)n, self->bases);
    self->bases[n] = '\0';
    plus = self->bases;
  }
  long differences = junction_flank_differences(self->refs, plus, (size_t)n,
      inner->ref_index, anchor, backwards);
  return (double)differences <= CHIMERA_MAX_DIFFERENCE * (double)n;
}

/*
   Whether the read stops at the match <m> (see chimera.h): the bases on
   the amplicon side of the match (a whole probe) are located in the
   references, and those on the other side neither continue the
   reference across the match (continues_reference()) nor are found next
   to the probe. The other side is given CHIMERA_OUTER_PROBES chances,
   the probe next to the match and those after it, as the segment can
   differ there from the references (e.g. at the cleavage site of the
   HA, just before its Uni12-like site): a read which runs through a
   site of the segment must not be cut there.
*/
static bool is_discontinuous(chimera_finder_t *self, const char *seq,
    size_t len, const terminus_match_t *m) {
  long probe = (long)reference_set_probe_len(self->refs);
  long before = m->from, after = m->to + 1;
  /* (the first base of the probe before the match) */
  long first = before - probe > 0 ? before - probe : 0;
  reference_hit_t inner, outer;
  /* (the amplicon side must hold a whole probe: fewer bases, e.g. those
     between the primer of the start of the read and a match just after
     it, are located too easily) */
  if (m->start ? after + probe > (long)len : before < probe)
    return false;
  if (m->start)
    locate(self, seq, len, after, after + probe, &inner);
  else
    locate(self, seq, len, before - probe, before, &inner);
  if (!inner.mapped)
    return false;
  if (continues_reference(self, seq, len, m, &inner))
    return false;
  for (long i = 0; i < CHIMERA_OUTER_PROBES; i++) {
    /* (the probe of the other side, <i> probes away from the match) */
    long tolerance = (i + 1) * CHIMERA_CONTINUITY_TOLERANCE;
    bool continuous;
    if (m->start) {
      long from = before - (i + 1) * probe;
      locate(self, seq, len, from, from + probe, &outer);
      continuous = is_continuous(&outer, &inner,
                                 after - (from > 0 ? from : 0), tolerance);
    } else {
      long from = after + i * probe;
      locate(self, seq, len, from, from + probe, &outer);
      continuous = is_continuous(&inner, &outer, from - first, tolerance);
    }
    if (continuous)
      return false;
  }
  return true;
}

/* whether the end of an amplicon <end> and the start of one <start>
   after it can be a pair in a read of <len> bases (see chimera.h) */
static bool is_pair(const terminus_match_t *end,
    const terminus_match_t *start, size_t len) {
  return !end->start && start->start &&
         start->from - end->to - 1 <= CHIMERA_MAX_GAP &&
         end->n_errors + start->n_errors <= CHIMERA_MAX_PAIR_ERRORS &&
         /* (bases enough for an amplicon on each side) */
         end->from >= CHIMERA_MIN_PART &&
         (long)len - 1 - start->to >= CHIMERA_MIN_PART;
}

/* whether the read runs through the pair <end>, <start>: the bases
   before it and those after it are next to each other in a reference */
static bool runs_through(chimera_finder_t *self, const char *seq,
    size_t len, const terminus_match_t *end, const terminus_match_t *start) {
  long probe = (long)reference_set_probe_len(self->refs);
  long from = end->from - probe > 0 ? end->from - probe : 0;
  reference_hit_t a, b;
  locate(self, seq, len, from, end->from, &a);
  locate(self, seq, len, start->to + 1, start->to + 1 + probe, &b);
  return is_continuous(&a, &b, start->to + 1 - from,
                       CHIMERA_OUTER_PROBES * CHIMERA_CONTINUITY_TOLERANCE);
}

/* end the current part at <to> and begin the next one at <from>; false
   if there is no room for another part */
static bool cut(chimera_part_t *parts, size_t *n, long to, long from) {
  if (*n + 1 >= CHIMERA_MAX_PARTS)
    return false;
  chimera_part_t *last = parts + *n - 1;
  /* (a part without any base, e.g. before a match at the start of the
     read, is not kept) */
  if (to < last->from)
    last->from = from;
  else {
    last->to = to;
    parts[(*n)++].from = from;
  }
  return true;
}

/* cut at the single match <m> if the read stops there: on the outer
   side of its primer */
static bool cut_single(chimera_finder_t *self, const char *seq, size_t len,
    const terminus_match_t *m, chimera_part_t *parts, size_t *n) {
  /* (not the primer of an end of the read) */
  bool inside = m->start ? m->from > self->window :
                m->to < (long)len - 1 - self->window;
  if (!inside || !is_discontinuous(self, seq, len, m))
    return true;
  long bound = terminus_matches_primer_bound(self->termini, seq, len, m);
  return m->start ? cut(parts, n, bound - 1, bound) :
         cut(parts, n, bound, bound + 1);
}

size_t chimera_finder_run(chimera_finder_t *self, const char *seq,
    size_t len, chimera_part_t *parts) {
  assert(self != NULL);
  assert(seq != NULL || len == 0);
  assert(parts != NULL);
  if (len == 0)
    return 0;
  terminus_match_t matches[TERMINUS_MATCHES_MAX];
  size_t n_matches = terminus_matches_find(self->termini, seq, len,
                                           matches);
  size_t n = 1;
  parts[0].from = 0;
  for (size_t i = 0; i < n_matches; i++) {
    const terminus_match_t *m = matches + i;
    bool ok = true;
    if (i + 1 < n_matches && is_pair(m, m + 1, len)) {
      /* a pair: the gap is dropped, unless the read runs through it
         (the two matches are then template) */
      if (self->refs == NULL || !runs_through(self, seq, len, m, m + 1))
        ok = cut(parts, &n, m->to, m[1].from);
      i++;
    } else if (self->refs != NULL) {
      ok = cut_single(self, seq, len, m, parts, &n);
    }
    if (!ok)
      break;
  }
  parts[n - 1].to = (long)len - 1;
  if (parts[n - 1].to < parts[n - 1].from)
    n--;
  return n;
}

/* --- unit tests --- */

/* split <read> and check the parts (<n> pairs of from and to) */
static void check_parts(chimera_finder_t *f, const char *read, size_t n,
    const long *expected) {
  chimera_part_t parts[CHIMERA_MAX_PARTS];
  test_eq_size(chimera_finder_run(f, read, strlen(read), parts), n);
  for (size_t i = 0; i < n; i++) {
    test_eq_long(parts[i].from, expected[2 * i]);
    test_eq_long(parts[i].to, expected[2 * i + 1]);
  }
}

void chimera_test(void) {
  Logger *logger = new_logger(stderr, LOG_ERROR);
  scheme_t *sc = scheme_new_builtin(SCHEME_DEFAULT, logger);
  enum { BODY = 200, SEG = 2 * BODY + 12 };
  const char *uni12 = "AGCAAAAGCAGG";
  const char *site = "AGCTATAGCAGG";      /* a site similar to Uni12 */
  const char *end13 = "CCTTGTTTCTACT";    /* rc of Uni13 */
  const char *start13 = "AGTAGAAACAAGG";  /* Uni13 */
  const char *gap = "ACGTTGCAACGTTGCAACGT";
  unsigned long state = 11;
  char a[BODY + 1], b[BODY + 1], other[BODY + 1], rc_a[BODY + 1];
  test_random_dna_seq(a, BODY, &state);
  test_random_dna_seq(b, BODY, &state);
  test_random_dna_seq(other, BODY, &state);
  primers_reverse_complement(a, BODY, rc_a);
  rc_a[BODY] = '\0';
  /* the segment: two halves around the site */
  char segment[SEG + 1];
  snprintf(segment, sizeof(segment), "%s%s%s", a, site, b);
  fasta_record_t rec = {"s", "", segment, SEG};
  fasta_t fa = {&rec, 1};
  reference_index_t *ix = reference_index_new(&fa);
  reference_set_t *rs = reference_set_new(ix, 32, 0.2);
  chimera_finder_t *f = chimera_finder_new(sc, 0.2, 100);
  char read[1024];
  long whole[2] = {0, 0};

  /* one amplicon: not split */
  snprintf(read, sizeof(read), "%s%s%s", uni12, a, end13);
  whole[1] = (long)strlen(read) - 1;
  check_parts(f, read, 1, whole);

  /* two amplicons: cut after the end of the first and before the start
     of the second, on either strand */
  snprintf(read, sizeof(read), "%s%s%s%s%s%s", uni12, a, end13, gap,
           uni12, b);
  long two[4] = {0, 12 + BODY + 13 - 1, 12 + BODY + 13 + 20,
                 (long)strlen(read) - 1};
  check_parts(f, read, 2, two);
  snprintf(read, sizeof(read), "%s%s%s%s%s%s", uni12, a, end13, gap,
           start13, b);
  two[3] = (long)strlen(read) - 1;
  check_parts(f, read, 2, two);

  /* a primer inside a read, without the end of an amplicon before it:
     a junction only with references */
  snprintf(read, sizeof(read), "%s%s%s%s", uni12, other, uni12, b);
  whole[1] = (long)strlen(read) - 1;
  check_parts(f, read, 1, whole);
  chimera_finder_use_references(f, rs);
  long single[4] = {0, 12 + BODY - 1, 12 + BODY, (long)strlen(read) - 1};
  check_parts(f, read, 2, single);

  /* the end of an amplicon on the other strand, then other bases */
  snprintf(read, sizeof(read), "%s%s%s", rc_a, "CCTGCTTTTGCT", other);
  long after[4] = {0, BODY + 12 - 1, BODY + 12, (long)strlen(read) - 1};
  check_parts(f, read, 2, after);

  /* the read runs through the site of the segment: not split */
  snprintf(read, sizeof(read), "%s%s", uni12, segment);
  whole[1] = (long)strlen(read) - 1;
  check_parts(f, read, 1, whole);

  /* the same, the bases before the site differing from the reference:
     those further away still tell that the read runs through it */
  memcpy(read + 12 + BODY - 24, "ACGTACGTACGTACGTACGTACGT", 24);
  check_parts(f, read, 1, whole);

  /* a read of a complete amplicon which runs through the site, of a
     sample far from the reference: a probe of the bases before the site
     is not located (one base in four differs), while the bases still
     fit the reference there, and the read is not split */
  snprintf(read, sizeof(read), "%s%s", uni12, segment);
  whole[1] = (long)strlen(read) - 1;
  for (long i = 12; i < 12 + BODY; i += 4)
    read[i] = read[i] == 'A' ? 'C' : 'A';
  check_parts(f, read, 1, whole);

  /* a read with an internal deletion whose junction holds the site:
     the bases before it lie farther from those after it on the
     reference, in the same order, so the read is one template */
  snprintf(read, sizeof(read), "%s%.100s%s%s", uni12, a, site, b + 100);
  whole[1] = (long)strlen(read) - 1;
  check_parts(f, read, 1, whole);

  /* a pair at the start of a read: the bases before it are no amplicon
     (the primer of the read end is one of the two matches) */
  snprintf(read, sizeof(read), "%.8s%s%s%s", other, end13, uni12, b);
  whole[1] = (long)strlen(read) - 1;
  check_parts(f, read, 1, whole);

  /* the start of an amplicon with a tailed primer (the reverse
     complement of the right end): cut before the tail, which goes with
     its primer */
  char rc_b[BODY + 1];
  primers_reverse_complement(b, BODY, rc_b);
  rc_b[BODY] = '\0';
  snprintf(read, sizeof(read), "%s%s%s%s", uni12, other,
           "ACGCGTGATCAGTAGAAACAAGG", rc_b);
  long tailed[4] = {0, 12 + BODY - 1, 12 + BODY, (long)strlen(read) - 1};
  check_parts(f, read, 2, tailed);

  /* a match just after the primer of the start of the read: the few
     bases before it do not tell where the read comes from */
  snprintf(read, sizeof(read), "%s%.8s%s%s", uni12, a, "CCTGCTTTTGCT",
           other);
  whole[1] = (long)strlen(read) - 1;
  check_parts(f, read, 1, whole);

  /* a pair with errors: two matches next to each other by chance are
     common with more than CHIMERA_MAX_PAIR_ERRORS errors together */
  chimera_finder_t *loose = chimera_finder_new(sc, 0.3, 100);
  snprintf(read, sizeof(read), "%s%s%s%s%s%s", uni12, other,
           "CCTTCTTTCTACT", gap, "AGCAAATGCAGG", other);
  long noisy[4] = {0, 12 + BODY + 13 - 1, 12 + BODY + 13 + 20,
                   (long)strlen(read) - 1};
  check_parts(loose, read, 2, noisy);
  snprintf(read, sizeof(read), "%s%s%s%s%s%s", uni12, other,
           "CCTTCTTACTACT", gap, "AGCTAATGCAGG", other);
  whole[1] = (long)strlen(read) - 1;
  check_parts(loose, read, 1, whole);
  chimera_finder_delete(loose);

  /* a pair in a read which runs through it (a segment which holds the
     two matches): not split with references, which tell */
  enum { SEG2 = 2 * BODY + 13 + 20 + 12 };
  char segment2[SEG2 + 1];
  snprintf(segment2, sizeof(segment2), "%s%s%s%s%s", a, end13, gap, uni12,
           b);
  fasta_record_t rec2 = {"t", "", segment2, SEG2};
  fasta_t fa2 = {&rec2, 1};
  reference_index_t *ix2 = reference_index_new(&fa2);
  reference_set_t *rs2 = reference_set_new(ix2, 32, 0.2);
  chimera_finder_t *f2 = chimera_finder_new(sc, 0.2, 100);
  snprintf(read, sizeof(read), "%s%s", uni12, segment2);
  long through[4] = {0, 12 + BODY + 13 - 1, 12 + BODY + 13 + 20,
                     (long)strlen(read) - 1};
  check_parts(f2, read, 2, through);
  chimera_finder_use_references(f2, rs2);
  whole[1] = (long)strlen(read) - 1;
  check_parts(f2, read, 1, whole);
  chimera_finder_delete(f2);
  reference_set_delete(rs2);
  reference_index_delete(ix2);

  /* an empty read */
  chimera_part_t parts[CHIMERA_MAX_PARTS];
  test_eq_size(chimera_finder_run(f, "", 0, parts), 0);

  chimera_finder_delete(f);
  reference_set_delete(rs);
  reference_index_delete(ix);
  scheme_delete(sc);
  delete_logger(logger);
}
