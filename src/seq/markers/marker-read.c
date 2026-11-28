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

/* read one marker out of a read */
static void read_marker(const marker_reference_t *site, size_t marker,
    const char *read, size_t len, double max_error_rate,
    marker_counts_t *counts) {
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
    if (span < 0)
      return;
    marker_span_value(span, value);
    marker_counts_add(counts, marker, value);
    return;
  }
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
    const char *read, size_t len, double max_error_rate,
    marker_counts_t *counts) {
  assert(chosen != NULL);
  assert(counts != NULL);
  if (read == NULL || len == 0)
    return;
  for (size_t m = 0; m < MARKER_N_MARKERS; m++) {
    if (chosen[m] == NULL)
      continue;
    if (ref >= 0 && (size_t)ref != chosen[m]->anchor.ref)
      continue;
    read_marker(chosen[m], m, read, len, max_error_rate, counts);
  }
}
