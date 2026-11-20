/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <limits.h>
#include <string.h>
#include "read-cover.h"
#include "reference.h"
#include "segments.h"
#include "test.h"

/* the segment (1 to 8, 0: unknown) of the reference <i> */
static int ref_segment(const fasta_t *refs, size_t i) {
  const char *description = refs->records[i].description;
  return description == NULL ? 0 : segments_of_description(description);
}

/* the distance of the length of the reference <i> from <span> */
static long length_distance(const fasta_t *refs, size_t i, long span) {
  long d = (long)refs->records[i].seq_len - span;
  return d < 0 ? -d : d;
}

/* the reference whose length is nearest to <span> (a complete amplicon's
   template length), among those considered by <set> (NULL: all); -1 if
   none is within LENGTH_TOLERANCE, or if a reference of another segment
   is too close in length to tell the two apart (e.g. PB2 and PB1; the
   references of a segment, e.g. of several subtypes, have similar
   lengths: the nearest is taken). A reference of unknown segment is
   taken as a segment of its own */
#define LENGTH_TOLERANCE 80
static long segment_by_length(const fasta_t *refs,
    const reference_set_t *set, long span) {
  long best = -1, best_d = LONG_MAX;
  for (size_t i = 0; i < refs->n_records; i++) {
    if (set != NULL && !reference_set_allows(set, i))
      continue;
    long d = length_distance(refs, i, span);
    if (d < best_d) {
      best_d = d;
      best = (long)i;
    }
  }
  if (best < 0 || best_d > LENGTH_TOLERANCE)
    return -1;
  int segment = ref_segment(refs, (size_t)best);
  for (size_t i = 0; i < refs->n_records; i++) {
    if (i == (size_t)best || (set != NULL && !reference_set_allows(set, i)))
      continue;
    if (segment != 0 && ref_segment(refs, i) == segment)
      continue;
    if (length_distance(refs, i, span) - best_d < LENGTH_TOLERANCE)
      return -1;   /* ambiguous */
  }
  return best;
}

/* the distance of the 5' end of the amplicon of a placed read from the
   terminus of the reference at which it would start; <minus>: the read
   is the minus strand of its amplicon (the amplicon then runs along the
   reference in the direction opposite to that of the read) */
static long placed_start_distance(const placement_t *p, bool minus) {
  return p->reverse != minus ? p->ref_len - p->to : p->from;
}

/* the class of a read with a 5' primer: terminal or internal, on its
   strand */
static size_t primed_cover(bool internal, bool minus) {
  if (internal)
    return minus ? COVER_INTERNAL_MINUS : COVER_INTERNAL_PLUS;
  return minus ? COVER_TERMINAL_MINUS : COVER_TERMINAL_PLUS;
}

/* the window within which the 5' end of a read of <span> bases, inferred
   from another position, is at the terminus of its segment: it grows
   with the span to allow for the indels of a long read */
static long start_window(long span) {
  return 60 + span / 20;
}

/* the span and the class of a read by the position <lp> of its left
   (5') end, and that of its right end <rp> if on the same reference
   (<both>). A position reported for an end read towards the start of the
   reference (a reverse hit, see reference_hit_t) is that of the last
   base of the span */
static void locate_by_left(const amplicon_t *a, const reference_hit_t *lp,
    const reference_hit_t *rp, bool both, long span, read_cover_t *c) {
  const amplicon_end_t *l = &a->left;
  bool minus = a->strand == STRAND_MINUS;
  c->ref = lp->ref_index;
  c->from = lp->ref_start;
  if (both) {
    /* the two ends: first and last base, in either order */
    c->to = rp->ref_start;
    if (c->from > c->to) {
      long t = c->from; c->from = c->to; c->to = t;
    }
    c->to++;
  } else if (lp->reverse) {
    /* the amplicon is on the minus strand of the reference: it runs
       from its 5' end towards the start of the reference */
    c->to = c->from + 1;
    c->from = c->to - span;
  } else {
    c->to = c->from + span;
  }
  /* the origin of the 5' end: that of its probe; without a primer the
     read is a fragment of unknown origin (grey); with a primer whose
     probe did not map, from where the end found from inside the read
     falls (the window allows for its approximation) */
  if (!l->found)
    c->cover = COVER_NO_PRIMER;
  else
    c->cover = primed_cover(l->ref.mapped ? l->origin == ORIGIN_INTERNAL :
        !reference_hit_is_terminal(lp, start_window(span)), minus);
}

/* the span and the class of a read of which only the 3' end is located,
   at <rp> (the 5' end was not found from inside the read either): the
   read is placed by it and its length, towards the start of the
   reference from the 3' end or, for an amplicon on the minus strand of
   the reference (the probe of its 3' end maps forward), towards the end
   of the reference. With no 5' primer it is a fragment of unknown origin
   (grey); with a 5' primer whose probe did not map (a noisy read start),
   the 5' origin is inferred from where the 5' end falls: within a window
   of the segment terminus is a full-length (terminal) amplicon, deeper
   is internal priming */
static void locate_by_right(const amplicon_t *a, const reference_hit_t *rp,
    long span, read_cover_t *c) {
  c->ref = rp->ref_index;
  if (rp->reverse) {
    c->to = rp->ref_start + 1;
    c->from = c->to - span;
  } else {
    c->from = rp->ref_start;
    c->to = c->from + span;
  }
  if (!a->left.found) {
    c->cover = COVER_NO_PRIMER;
    return;
  }
  long start = rp->reverse ? c->from : rp->ref_len - c->to;
  c->cover = primed_cover(start > start_window(span),
                          a->strand == STRAND_MINUS);
}

/* the span and the class of a read of which no end mapped, placed by a
   probe from inside it (see placement.h); its 5' origin as above */
static void locate_by_placement(const amplicon_t *a, long span,
    read_cover_t *c) {
  const placement_t *p = &a->placement;
  bool minus = a->strand == STRAND_MINUS;
  c->ref = p->ref_index;
  c->from = p->from;
  c->to = p->to;
  c->cover = !a->left.found ? COVER_NO_PRIMER :
      primed_cover(placed_start_distance(p, minus) > start_window(span),
                   minus);
}

/* a read with both ends mapped much farther apart than its template:
   its two flanks, up to its junction if located, or else half of the
   template on each side */
static void set_deletion(const amplicon_t *a, long span, read_cover_t *c) {
  c->cover = COVER_DELETION;
  if (a->junction.found) {
    c->left_to = a->junction.left_to;
    c->right_from = a->junction.right_from;
  } else {
    c->left_to = c->from + span / 2;
    c->right_from = c->to - (span - span / 2);
  }
}

bool read_cover_locate(const fasta_t *refs, const reference_set_t *set,
    const amplicon_t *a, read_cover_t *cover) {
  assert(a != NULL);
  assert(cover != NULL);
  memset(cover, 0, sizeof(read_cover_t));
  if (refs == NULL)
    return false;
  const amplicon_end_t *l = &a->left, *r = &a->right;
  /* the positions of the ends: those of their probes, or else those
     found from inside the read (see amplicon.h) */
  const reference_hit_t *lp = amplicon_end_position(l),
                        *rp = amplicon_end_position(r);
  long span = (long)amplicon_kept_length(a), ref;
  bool both = false;
  if (lp != NULL) {
    both = rp != NULL && rp->ref_index == lp->ref_index;
    locate_by_left(a, lp, rp, both, span, cover);
  } else if (rp != NULL) {
    locate_by_right(a, rp, span, cover);
  } else if (a->placement.found) {
    locate_by_placement(a, span, cover);
  } else if (l->found && r->found && span > 0 &&
             (ref = segment_by_length(refs, set, span)) >= 0) {
    /* a complete amplicon whose probe did not map (e.g. a noisy read
       start): placed over the whole segment of its length, as terminal */
    cover->ref = (size_t)ref;
    cover->from = 0;
    cover->to = (long)refs->records[ref].seq_len;
    cover->cover = primed_cover(false, a->strand == STRAND_MINUS);
  } else {
    return false;   /* not located and not placeable by length */
  }
  if (both && (a->junction.found ||
               junction_expected(span, cover->to - cover->from)))
    set_deletion(a, span, cover);
  cover->located = true;
  return true;
}

size_t read_cover_category(const read_cover_t *cover) {
  assert(cover != NULL);
  if (!cover->located)
    return READ_CATEGORY_OTHER;
  switch (cover->cover) {
    case COVER_TERMINAL_PLUS:
    case COVER_TERMINAL_MINUS:
      return READ_CATEGORY_FULL;
    case COVER_INTERNAL_PLUS:
    case COVER_INTERNAL_MINUS:
      return READ_CATEGORY_INTERNAL;
    case COVER_DELETION:
      return READ_CATEGORY_DELETION;
    default:
      return READ_CATEGORY_OTHER;
  }
}

const char* read_category_name(size_t category) {
  static const char *const names[READ_N_CATEGORIES] = {
      "full", "internal", "deletion", "other"};
  assert(category < READ_N_CATEGORIES);
  return names[category];
}

/* --- unit tests --- */

/* an amplicon of <span> bases with both primers found, no end located */
static amplicon_t test_amplicon(long span) {
  amplicon_t a;
  memset(&a, 0, sizeof(a));
  a.strand = STRAND_PLUS;
  a.left.found = a.right.found = true;
  a.keep_from = 30;
  a.keep_to = 30 + span - 1;
  return a;
}

/* the end <end> mapped at <pos> of the reference 0 (of 2000 bases) */
static void test_map(amplicon_end_t *end, long pos, end_origin_t origin) {
  end->origin = origin;
  end->ref.mapped = true;
  end->ref.ref_index = 0;
  end->ref.ref_len = 2000;
  end->ref.ref_start = pos;
}

/* the category of <a> on <refs> */
static size_t test_category(const fasta_t *refs, const amplicon_t *a,
    read_cover_t *c) {
  read_cover_locate(refs, NULL, a, c);
  return read_cover_category(c);
}

void read_cover_test(void) {
  static char seq[2001], seq2[1001];
  memset(seq, 'A', 2000);
  memset(seq2, 'A', 1000);
  fasta_record_t recs[2] = {{"r0", "", seq, 2000}, {"r1", "", seq2, 1000}};
  fasta_t refs = {recs, 2};
  read_cover_t c;
  test_eq_str(read_category_name(READ_CATEGORY_FULL), "full");
  test_eq_str(read_category_name(READ_CATEGORY_OTHER), "other");
  /* no references, or nothing located: other */
  amplicon_t a = test_amplicon(300);
  test_assert(!read_cover_locate(NULL, NULL, &a, &c));
  test_eq_size(read_cover_category(&c), READ_CATEGORY_OTHER);
  test_assert(!read_cover_locate(&refs, NULL, &a, &c));
  test_eq_size(read_cover_category(&c), READ_CATEGORY_OTHER);
  /* a complete amplicon of the length of a segment: a full template */
  a = test_amplicon(1050);
  test_eq_size(test_category(&refs, &a, &c), READ_CATEGORY_FULL);
  test_eq_size(c.ref, 1);
  test_eq_long(c.from, 0);
  test_eq_long(c.to, 1000);
  /* not if its length is too far from that of any segment, or as near
     to that of two */
  a = test_amplicon(1100);
  test_eq_size(test_category(&refs, &a, &c), READ_CATEGORY_OTHER);
  fasta_record_t close[2] = {{"r0", "", seq, 1040}, {"r1", "", seq2, 1000}};
  fasta_t close_refs = {close, 2};
  a = test_amplicon(1010);
  test_eq_size(test_category(&close_refs, &a, &c), READ_CATEGORY_OTHER);
  /* unless they are of the same segment: the nearest */
  close[0].description = "A/x/1/2000(H1N1) segment 4";
  close[1].description = "A/y/2/2001(H3N2) segment 4";
  test_eq_size(test_category(&close_refs, &a, &c), READ_CATEGORY_FULL);
  test_eq_size(c.ref, 1);
  close[1].description = "A/y/2/2001(H3N2) segment 6";
  test_eq_size(test_category(&close_refs, &a, &c), READ_CATEGORY_OTHER);
  /* both ends mapped, as far apart as the template: full, by strand */
  a = test_amplicon(2000);
  test_map(&a.left, 0, ORIGIN_TERMINAL);
  test_map(&a.right, 1999, ORIGIN_TERMINAL);
  test_eq_size(test_category(&refs, &a, &c), READ_CATEGORY_FULL);
  test_eq_size(c.cover, COVER_TERMINAL_PLUS);
  a.strand = STRAND_MINUS;
  test_eq_size(test_category(&refs, &a, &c), READ_CATEGORY_FULL);
  test_eq_size(c.cover, COVER_TERMINAL_MINUS);
  /* much farther apart than the template: a deletion, half of the
     template on each side, or the flanks up to the junction */
  a = test_amplicon(400);
  test_map(&a.left, 0, ORIGIN_TERMINAL);
  test_map(&a.right, 1999, ORIGIN_TERMINAL);
  test_eq_size(test_category(&refs, &a, &c), READ_CATEGORY_DELETION);
  test_eq_long(c.from, 0);
  test_eq_long(c.left_to, 200);
  test_eq_long(c.right_from, 1800);
  test_eq_long(c.to, 2000);
  a.junction = (junction_t){.found = true, .left_to = 100,
                            .right_from = 1700};
  test_eq_size(test_category(&refs, &a, &c), READ_CATEGORY_DELETION);
  test_eq_long(c.left_to, 100);
  test_eq_long(c.right_from, 1700);
  /* the left end internal: internal priming, forward and reverse */
  a = test_amplicon(300);
  test_map(&a.left, 1000, ORIGIN_INTERNAL);
  test_eq_size(test_category(&refs, &a, &c), READ_CATEGORY_INTERNAL);
  test_eq_size(c.cover, COVER_INTERNAL_PLUS);
  test_eq_long(c.from, 1000);
  test_eq_long(c.to, 1300);
  a.left.ref.reverse = true;
  test_eq_size(test_category(&refs, &a, &c), READ_CATEGORY_INTERNAL);
  test_eq_long(c.from, 701);
  test_eq_long(c.to, 1001);
  /* located by an end, but without a 5' primer: other */
  a.left.found = false;
  a.left.inner = a.left.ref;
  test_eq_size(test_category(&refs, &a, &c), READ_CATEGORY_OTHER);
  test_assert(c.located);
  test_eq_size(c.cover, COVER_NO_PRIMER);
  /* the right end only: the 5' end is where the length puts it, here
     far from the end of the reference (internal), then at it (full) */
  a = test_amplicon(300);
  test_map(&a.right, 100, ORIGIN_TERMINAL);
  test_eq_size(test_category(&refs, &a, &c), READ_CATEGORY_INTERNAL);
  test_eq_long(c.from, 100);
  test_eq_long(c.to, 400);
  a.right.ref.ref_start = 1690;
  test_eq_size(test_category(&refs, &a, &c), READ_CATEGORY_FULL);
  a.left.found = false;
  test_eq_size(test_category(&refs, &a, &c), READ_CATEGORY_OTHER);
  /* placed by a probe from inside the read */
  a = test_amplicon(300);
  a.placement = (placement_t){.found = true, .ref_index = 0,
                              .ref_len = 2000, .from = 500, .to = 800};
  test_eq_size(test_category(&refs, &a, &c), READ_CATEGORY_INTERNAL);
  a.placement.from = 10;
  a.placement.to = 310;
  test_eq_size(test_category(&refs, &a, &c), READ_CATEGORY_FULL);
  a.left.found = false;
  test_eq_size(test_category(&refs, &a, &c), READ_CATEGORY_OTHER);
}
