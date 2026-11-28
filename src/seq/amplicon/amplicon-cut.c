/* (c) 2026 Giorgio Gonnella */

#include <stdbool.h>
#include <stddef.h>
#include "amplicon-analyzer.h"
#include "amplicon-cut.h"
#include "primers.h"
#include "reference.h"

/* a value limited to [lo, hi] */
static long clamp(long v, long lo, long hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}

/* the segment of the reference the hit <hit> maps to (0: unknown) */
static int hit_segment(const amplicon_analyzer_t *self,
    const reference_hit_t *hit) {
  if (!hit->mapped || self->ref_segment == NULL)
    return 0;
  return self->ref_segment[hit->ref_index];
}

/*
   Locate the probe <probe> (length <len>) of an end in the references,
   and set the origin of the end accordingly; <primer> is the primer
   found at that end, whose terminus is kept only at a terminal end of
   a segment it binds (see primers.h): a primer which reaches past the
   universal terminus into segment-specific bases writes its own
   sequence at the end of any other segment, as at an internally primed
   end, so those bases are removed there.
*/
static void classify_end(amplicon_analyzer_t *self, const primer_t *primer,
    const char *probe, size_t len, amplicon_end_t *end) {
  reference_set_locate(self->refs, probe, len, &end->ref);
  if (!end->ref.mapped)
    end->origin = ORIGIN_UNKNOWN;
  else
    end->origin = reference_hit_is_terminal(&end->ref,
        self->terminal_window) ? ORIGIN_TERMINAL : ORIGIN_INTERNAL;
  end->keep_terminus = end->origin == ORIGIN_TERMINAL &&
      primers_binds_segment(primer, hit_segment(self, &end->ref));
}

/*
   Classify the ends of the amplicon <a> found in <seq> (length <len>)
   against the references (see amplicon.h): the probe of the left end is
   the sequence after its primer; that of the right end the reverse
   complement of the sequence before its primer, so that both are read
   from the primer into the template, as the strand the primer started.
   The probes are taken from the bases before the other primer only.
*/
void amplicon_classify_ends(amplicon_analyzer_t *self, const char *seq,
    size_t len, amplicon_t *a) {
  long lo = a->left.found ? a->left.hit.end + 1 : 0;
  long hi = a->right.found ? a->right.hit.start : (long)len;
  long avail = hi > lo ? hi - lo : 0;
  long probe_len = (long)reference_set_probe_len(self->refs);
  long n = avail < probe_len ? avail : probe_len;
  if (a->left.found)
    classify_end(self, self->scheme->left + a->left.variant, seq + lo,
                 (size_t)n, &a->left);
  if (a->right.found) {
    primers_reverse_complement(seq + hi - n, (size_t)n, self->probe);
    self->probe[n] = '\0';
    classify_end(self, self->scheme->right + a->right.variant,
                 self->probe, (size_t)n, &a->right);
  }
}

/*
   The part of <seq> (length <len>) to keep, in the orientation in which
   the layout was found; an end started by internal priming is trimmed
   as by TRIM_PRIMER, whatever the policy (see amplicon.h)
*/
void amplicon_compute_keep(const amplicon_t *a, trim_policy_t policy,
    size_t len, long *from, long *to) {
  *from = 0;
  *to = (long)len - 1;
  /* the conserved terminus is kept only when it is confidently terminal
     (TRIM_TERMINUS and ORIGIN_TERMINAL); an internal or an unlocated end
     ("cut when unsure") has its whole primer removed, so that a primer
     base never reaches the consensus unless the end is known to be a
     segment terminus */
  if (a->left.found)
    *from = policy == TRIM_TERMINUS && a->left.keep_terminus ?
            clamp(a->left.hit.terminus_start, 0, (long)len) :
            a->left.hit.end + 1;
  if (a->right.found)
    *to = policy == TRIM_TERMINUS && a->right.keep_terminus ?
          clamp(a->right.hit.terminus_end, 0, (long)len) - 1 :
          a->right.hit.start - 1;
  *from = clamp(*from, 0, (long)len);
  *to = clamp(*to, -1, (long)len - 1);
  /* no template between the two primers (a primer dimer): nothing is
     kept, whatever the policy (the termini alone are primer sequence) */
  if (a->left.found && a->right.found &&
      a->right.hit.start <= a->left.hit.end + 1) {
    *from = 1;
    *to = 0;
  }
}

/*
   A read can have a primer of the same end of the scheme at both of its
   ends, one of them reverse-complemented: two amplicons joined without
   a primer between them (see chimera.h), each read from its primer
   towards the junction. Only one of the two is an end of the layout
   found; the other, an end of the layout <other> of the opposite
   orientation, is removed whole with what lies outside it, as it is
   primer sequence all the same: <from> and <to>, the part of the read
   of length <len> to keep in the orientation of <a>, are limited
   accordingly.
*/
void amplicon_cut_facing_primer(const amplicon_t *a, const amplicon_t *other,
    size_t len, long *from, long *to) {
  long n = (long)len;
  if (!a->right.found && other->left.found) {
    /* (the positions of <other> are mirrored) */
    long start = n - 1 - other->left.hit.end;
    if ((!a->left.found || start > a->left.hit.end) && start - 1 < *to)
      *to = start - 1;
  }
  if (!a->left.found && other->right.found) {
    long end = n - 1 - other->right.hit.start;
    if ((!a->right.found || end < a->right.hit.start) && end + 1 > *from)
      *from = end + 1;
  }
}
