/* (c) 2026 Giorgio Gonnella */

#include <stdbool.h>
#include <stddef.h>
#include "alloc.h"
#include "amplicon-analyzer.h"
#include "amplicon-locate.h"
#include "junction.h"
#include "placement.h"
#include "reference.h"
#include "sequtil.h"

/* place on the references the kept part of a read of which no end was
   located (see "Placement" in amplicon.h) */
void amplicon_place(amplicon_analyzer_t *self, const char *seq,
    amplicon_t *a) {
  a->placement.found = false;
  if (self->refs == NULL || (a->left.found && a->left.ref.mapped) ||
      (a->right.found && a->right.ref.mapped))
    return;
  size_t kept = amplicon_kept_length(a);
  if (kept > 0)
    placement_locate(self->refs, seq + a->keep_from, kept, &a->placement);
}

/* the template of the read, in the analyzed orientation: the bases
   between the primers, or, at an end without a primer, from the end of
   the kept part; empty where nothing is kept (keep_from and keep_to
   are then a sentinel pair, not real positions of the read) */
static void template_span(const amplicon_t *a, size_t len, long *from,
    long *to) {
  if (amplicon_kept_length(a) == 0) {
    *from = 0;
    *to = -1;
    return;
  }
  bool minus = a->strand == STRAND_MINUS;
  *from = a->left.found ? a->left.hit.end + 1 :
          (minus ? (long)len - 1 - a->keep_to : a->keep_from);
  *to = a->right.found ? a->right.hit.start - 1 :
        (minus ? (long)len - 1 - a->keep_from : a->keep_to);
}

/* the greatest share of the bases between an end and the probe which
   located it from inside which may differ from the reference, for the
   bases to be taken as continuing the reference there, and the same at
   a terminus (where a sample can differ more from a reference) */
#define INNER_MAX_DIFFERENCE 0.3
#define INNER_MAX_TERMINAL_DIFFERENCE 0.4

/* where the template of a complete amplicon begins or ends in a
   reference: after its conserved terminus (at most this many bases),
   or, in a reference without its terminus, at its first or last base;
   a reference can lack some more bases, or have a few more: the end of
   a read is searched up to INNER_TERMINUS_SLACK bases from both */
#define INNER_TERMINUS_LEN 13
#define INNER_TERMINUS_SLACK 8

/*
   The position of an end with a primer, located by a probe <shift> bases
   inside the template, at <hit>: the <shift> bases <skipped> between the
   end and the probe (read from the end inwards) are compared with the
   reference. If they continue it up to the probe, the end is <shift>
   bases before the probe. If they do not, the template has a junction
   among them: the probe is of the other flank of an internal deletion,
   and the end is where the primer of a complete amplicon is, at the
   terminus of the segment, if the bases agree with the reference there;
   otherwise the end is taken as if the bases continued the reference.
*/
static long inner_end_position(amplicon_analyzer_t *self,
    const char *skipped, long shift, const reference_hit_t *hit) {
  long continued = hit->ref_start + (hit->reverse ? shift : -shift);
  if (shift == 0)
    return continued;
  /* the skipped bases on the plus strand of the reference: before the
     probe, or, for a probe of the minus strand, after it */
  const char *plus = skipped;
  if (hit->reverse) {
    XCHECKSPACEINC(self->template, &self->template_cap, (size_t)shift + 1,
                   1024);
    reverse_complement(skipped, (size_t)shift, self->template);
    self->template[shift] = '\0';
    plus = self->template;
  }
  long differences = junction_flank_differences(self->refs, plus,
      (size_t)shift, hit->ref_index,
      hit->reverse ? hit->ref_start + 1 : hit->ref_start - 1,
      !hit->reverse);
  if ((double)differences <= INNER_MAX_DIFFERENCE * (double)shift)
    return continued;
  long best = -1, best_end = continued;
  for (long d = -INNER_TERMINUS_SLACK;
       d <= INNER_TERMINUS_LEN + INNER_TERMINUS_SLACK; d++) {
    /* (d: bases of the reference before the end, at that terminus) */
    long end = hit->reverse ? hit->ref_len - 1 - d : d;
    differences = junction_flank_differences(self->refs, plus,
        (size_t)shift, hit->ref_index, end, hit->reverse);
    if (best < 0 || differences < best) {
      best = differences;
      best_end = end;
    }
  }
  return (double)best <= INNER_MAX_TERMINAL_DIFFERENCE * (double)shift ?
         best_end : continued;
}

/* locate from inside the template (<n> bases at <template> in the read)
   one of its ends, the one at its start (<at_start>) or the one at its
   end, by probes read from that end inwards; <other>: the position of
   the other end, if known: the probe must then map to its reference, on
   the other strand; <primer>: the end has a primer (see
   inner_end_position()). Into <out> (mapped: false if no probe is
   located) */
static void locate_from_inside(amplicon_analyzer_t *self,
    const char *template, long n, bool at_start, bool primer,
    const reference_hit_t *other, reference_hit_t *out) {
  long probe_len = (long)reference_set_probe_len(self->refs);
  out->mapped = false;
  for (long i = 0; i < AMPLICON_INNER_PROBES; i++) {
    /* (the probe must not reach the bases of the probe of the other
       end) */
    long shift = i * probe_len;
    if (shift + 2 * probe_len > n)
      break;
    const char *probe = template + shift;
    if (!at_start) {
      reverse_complement(template + n - shift - probe_len,
                         (size_t)probe_len, self->probe);
      self->probe[probe_len] = '\0';
      probe = self->probe;
    }
    /* with the other end: in its reference, on the other strand;
       without: wherever the k-mers of the probe find it */
    reference_hit_t hit;
    if (other != NULL)
      reference_set_locate_in(self->refs, probe, (size_t)probe_len,
                              other->ref_index, !other->reverse, &hit);
    else
      reference_set_locate_seeded(self->refs, probe, (size_t)probe_len,
                                  &hit);
    if (!hit.mapped)
      continue;
    /* from the first base of the probe out to the end of the template */
    if (!primer) {
      hit.ref_start += hit.reverse ? shift : -shift;
    } else {
      const char *skipped = template;
      if (!at_start && shift > 0) {
        XCHECKSPACEINC(self->skipped, &self->skipped_cap,
                       (size_t)shift + 1, 256);
        reverse_complement(template + n - shift, (size_t)shift,
                           self->skipped);
        self->skipped[shift] = '\0';
        skipped = self->skipped;
      }
      hit.ref_start = inner_end_position(self, skipped, shift, &hit);
    }
    *out = hit;
    return;
  }
}

/* locate from inside the read the ends which were not located (see
   "Ends located from inside" in amplicon.h): each with the position of
   the other one, if known; if neither is, the left one first, then the
   right one, then the left one again if only the right one was found */
void amplicon_locate_inner_ends(amplicon_analyzer_t *self, const char *seq,
    size_t len, amplicon_t *a) {
  a->left.inner.mapped = a->right.inner.mapped = false;
  if (self->refs == NULL)
    return;
  long from, to;
  template_span(a, len, &from, &to);
  long n = to - from + 1;
  if (n <= 0)
    return;
  bool minus = a->strand == STRAND_MINUS;
  /* the template in the read: the left end of the amplicon is its start,
     or, in a read of the minus strand, its end */
  const char *template = seq + (minus ? (long)len - 1 - to : from);
  bool left_known = amplicon_end_position(&a->left) != NULL,
       right_known = amplicon_end_position(&a->right) != NULL;
  if (amplicon_end_position(&a->left) == NULL)
    locate_from_inside(self, template, n, !minus, a->left.found,
                       amplicon_end_position(&a->right), &a->left.inner);
  if (!right_known)
    locate_from_inside(self, template, n, minus, a->right.found,
                       amplicon_end_position(&a->left), &a->right.inner);
  if (!right_known && amplicon_end_position(&a->left) == NULL &&
      amplicon_end_position(&a->right) != NULL)
    locate_from_inside(self, template, n, !minus, a->left.found,
                       amplicon_end_position(&a->right), &a->left.inner);
  /* a read of which no end was located by its own probe needs both
     ends: one probe alone can be a chance match of a read which is not
     of the references (such a read is left to its placement) */
  if (!left_known && !right_known &&
      (!a->left.inner.mapped || !a->right.inner.mapped))
    a->left.inner.mapped = a->right.inner.mapped = false;
}

/* locate the junction of a read whose two ends mapped to the same
   reference, farther apart than the template between them is long (see
   "Internal deletions" in amplicon.h) */
void amplicon_find_junction(amplicon_analyzer_t *self, const char *seq,
    size_t len, amplicon_t *a) {
  a->junction.found = false;
  if (self->refs == NULL)
    return;
  const reference_hit_t *l = amplicon_end_position(&a->left),
                        *r = amplicon_end_position(&a->right);
  if (l == NULL || r == NULL || l->ref_index != r->ref_index ||
      l->reverse == r->reverse)
    return;
  /* the template, in the analyzed orientation and in the read */
  long from, to;
  template_span(a, len, &from, &to);
  if (to < from)
    return;
  size_t n = (size_t)(to - from + 1);
  bool minus = a->strand == STRAND_MINUS;
  long template_from = minus ? (long)len - 1 - to : from;
  const char *template = seq + template_from;
  /* its first and last base on the plus strand of the reference */
  long first = l->reverse ? r->ref_start : l->ref_start;
  long last = l->reverse ? l->ref_start : r->ref_start;
  if (!junction_expected((long)n, last - first + 1))
    return;
  bool reversed = minus != l->reverse;
  if (reversed) {
    /* the read is the minus strand of the reference */
    XCHECKSPACEINC(self->template, &self->template_cap, n + 1, 1024);
    reverse_complement(template, n, self->template);
    self->template[n] = '\0';
    template = self->template;
  }
  junction_locate(self->refs, template, n, l->ref_index, first, last,
                  &a->junction);
  /* (in a read of the minus strand of the reference, the right flank
     comes first) */
  if (a->junction.found)
    a->junction_at = template_from + (reversed ?
        (long)n - a->junction.offset : a->junction.offset);
}
