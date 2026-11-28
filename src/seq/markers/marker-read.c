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

/* read one marker out of a read */
static void read_marker(const marker_reference_t *site, size_t marker,
    const char *read, size_t len, const marker_breaks_t *breaks,
    double max_error_rate, marker_counts_t *counts) {
  const marker_def_t *def = marker_def(marker);
  marker_counts_set_reference(counts, marker, site->value, site->anchor_aa);
  marker_counts_searched(counts, marker);
  anchor_hit_t hit;
  if (!anchor_find(&site->anchor, read, len, max_error_rate, &hit))
    return;
  /* the anchor as the read carries it: the reads which agree on it are
     the ones the marker beside it is read from */
  char nt[MARKER_WINDOW_NT + 1];
  char value[MARKER_VALUE_MAX + 1];
  if (anchor_bases_at(&hit, read, len, MARKER_ANCHOR_NT, nt)) {
    char aa[MARKER_ANCHOR_AA + 1];
    codons_translate(nt, MARKER_ANCHOR_NT, aa);
    marker_counts_add_anchor(counts, marker, aa);
  }
  if (def->kind == MARKER_SPAN) {
    anchor_hit_t hit2;
    if (!anchor_find(&site->anchor2, read, len, max_error_rate, &hit2))
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
