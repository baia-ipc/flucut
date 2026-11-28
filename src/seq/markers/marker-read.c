/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <string.h>
#include "codons.h"
#include "marker-read.h"

/*
   The bases between the first base of the anchor and the first base of
   the second one, counted in the orientation of the reference, or -1 if
   the two matches do not lie that way (the read matched one of them on
   the other strand, which does not say how far apart they are)
*/
static long span_of(const anchor_hit_t *first, const anchor_hit_t *second) {
  if (first->reverse != second->reverse)
    return -1;
  /* in a read which runs the other way, the first base of an anchor is
     the last base of its match, and the stretch runs backwards */
  long span = first->reverse ? first->end - second->end
                             : second->start - first->start;
  return span > 0 ? span : -1;
}

/* whether the read follows its reference from <from> to <to>: no cut
   falls inside that interval and no blocked stretch overlaps it */
static bool unbroken(const marker_breaks_t *breaks, long from, long to) {
  for (size_t i = 0; i < breaks->n_cuts; i++)
    if (breaks->cuts[i] > from && breaks->cuts[i] <= to)
      return false;
  for (size_t i = 0; i < breaks->n_blocked; i++)
    if (breaks->blocked[i].from <= to && breaks->blocked[i].to >= from)
      return false;
  return true;
}

/*
   Whether <at>, the position of the anchor in the reference, lies in a
   stretch the read covers. A read which was not located covers nothing
   known, and is accepted: its anchor is all that places it.
*/
static bool placed(const marker_breaks_t *breaks, long at) {
  if (breaks->n_covered == 0)
    return true;
  for (size_t i = 0; i < breaks->n_covered; i++)
    if (at >= breaks->covered[i].from - MARKER_PLACEMENT_SLACK &&
        at <= breaks->covered[i].to + MARKER_PLACEMENT_SLACK)
      return true;
  return false;
}

/*
   Find <anchor> in the read. A read located on its reference carries
   the anchor about where the stretches it covers put it: counted along
   them from the first base kept of the read if the read runs along the
   reference, from the last one if it runs the other way. The anchor is
   looked for within a margin of those two places first - the reads of a
   whole segment are long, and the anchor is short - and in the whole
   read only if it is not found there, so that a read whose place was
   estimated badly is still read. The margin allows for the error of the
   placement and for the indels of a long read.
*/
static bool find_anchor(const anchor_t *anchor, const char *read,
    size_t len, const marker_breaks_t *breaks, double max_error_rate,
    anchor_hit_t *hit) {
  long offset = -1, total = 0;
  for (size_t i = 0; i < breaks->n_covered; i++) {
    const marker_span_t *c = &breaks->covered[i];
    if (offset < 0 && anchor->ref_pos >= c->from - MARKER_PLACEMENT_SLACK &&
        anchor->ref_pos <= c->to + MARKER_PLACEMENT_SLACK)
      offset = total + (anchor->ref_pos - c->from);
    total += c->to - c->from + 1;
  }
  if (offset >= 0) {
    long alen = (long)anchor->len;
    long plus_at = breaks->read_from + offset;
    long minus_at = breaks->read_from + total - offset - alen;
    long plus_margin = MARKER_PLACEMENT_SLACK + offset / 10;
    long minus_margin = MARKER_PLACEMENT_SLACK + (total - offset) / 10;
    if (anchor_find_within(anchor, read, len, plus_at - plus_margin,
            plus_at + alen - 1 + plus_margin, minus_at - minus_margin,
            minus_at + alen - 1 + minus_margin, max_error_rate, hit))
      return true;
  }
  return anchor_find(anchor, read, len, max_error_rate, hit);
}

/* whether the anchor matched, and the window beside it was read, on a
   contiguous stretch of the read */
static bool window_whole(const anchor_hit_t *hit, size_t len,
    const marker_breaks_t *breaks, long offset, size_t n) {
  long from, to;
  if (!anchor_window_span(hit, len, offset, n, &from, &to))
    return false;
  if (hit->start < from)
    from = hit->start;
  if (hit->end > to)
    to = hit->end;
  return unbroken(breaks, from, to);
}

/*
   Whether the residues <read_aa> a read carries where the anchor
   matched and the residues <reference_aa> of the anchor in the
   reference have the same residue anywhere. The anchor is conserved, so
   a read which matched it carries almost all of it; a read which shares
   not one residue with it matched bases which translate to something
   else, because the match is in another frame or is not of the anchor.
   What lies beside such a match is not the position the marker names,
   so the read gives no value - the residues it carries at the anchor
   are still counted, because they are the evidence of it.

   Between two references of the same marker in the built-in set the
   fewest shared of these residues is two of eight (the
   receptor-binding anchor, between the most distant subtypes), so a
   strain difference does not reach none.
*/
static bool anchor_agrees(const char *read_aa, const char *reference_aa) {
  if (read_aa[0] == '\0' || reference_aa[0] == '\0')
    return true;
  for (size_t i = 0; read_aa[i] != '\0' && reference_aa[i] != '\0'; i++)
    if (read_aa[i] == reference_aa[i])
      return true;
  return false;
}

/* read one marker out of a read */
static void read_marker(const marker_reference_t *site, size_t marker,
    const char *read, size_t len, const marker_breaks_t *breaks,
    double max_error_rate, marker_counts_t *counts) {
  const marker_def_t *def = marker_def(marker);
  marker_counts_set_reference(counts, marker, site->value, site->anchor_aa);
  marker_counts_searched(counts, marker);
  if (!placed(breaks, site->anchor.ref_pos))
    return;
  anchor_hit_t hit;
  if (!find_anchor(&site->anchor, read, len, breaks, max_error_rate, &hit))
    return;
  /* the anchor as the read carries it: the reads which agree on it are
     the ones the marker beside it is read from */
  char nt[MARKER_WINDOW_NT + 1];
  char value[MARKER_VALUE_MAX + 1];
  if (anchor_bases_at(&hit, read, len, MARKER_ANCHOR_NT, nt)) {
    char aa[MARKER_ANCHOR_AA + 1];
    codons_translate(nt, MARKER_ANCHOR_NT, aa);
    marker_counts_add_anchor(counts, marker, aa);
    if (!anchor_agrees(aa, site->anchor_aa))
      return;
  }
  if (def->kind == MARKER_SPAN) {
    anchor_hit_t hit2;
    if (!placed(breaks, site->anchor2.ref_pos) ||
        !find_anchor(&site->anchor2, read, len, breaks, max_error_rate,
                     &hit2))
      return;
    long span = span_of(&hit, &hit2);
    if (span < 0 ||
        !unbroken(breaks,
                  hit.start < hit2.start ? hit.start : hit2.start,
                  hit.end > hit2.end ? hit.end : hit2.end))
      return;
    marker_span_value(span, value);
    marker_counts_add(counts, marker, value);
    return;
  }
  long offset = def->offset < 0 ? def->offset
                                : (long)anchor_hit_len(&hit) + def->offset;
  if (!window_whole(&hit, len, breaks, offset, def->window))
    return;
  bool got = def->offset < 0 ?
      anchor_bases_window(&hit, read, len, def->offset, def->window, nt) :
      anchor_bases_after(&hit, read, len, (size_t)def->offset, def->window,
                         nt);
  if (!got)
    return;
  if (marker_value_of(def->kind, nt, def->window, value))
    marker_counts_add(counts, marker, value);
}

void marker_read_sites(const marker_reference_t *const *chosen, long ref,
    const marker_breaks_t *breaks, const char *read, size_t len,
    double max_error_rate, marker_counts_t *counts) {
  assert(chosen != NULL);
  assert(breaks != NULL);
  assert(counts != NULL);
  if (read == NULL || len == 0)
    return;
  for (size_t m = 0; m < MARKER_N_MARKERS; m++) {
    if (chosen[m] == NULL)
      continue;
    if (ref >= 0 && (size_t)ref != chosen[m]->anchor.ref)
      continue;
    read_marker(chosen[m], m, read, len, breaks, max_error_rate, counts);
  }
}
