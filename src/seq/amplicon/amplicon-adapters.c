/* (c) 2026 Giorgio Gonnella */

#include <stdbool.h>
#include <stddef.h>
#include "adapters.h"
#include "amplicon-adapters.h"
#include "amplicon-analyzer.h"

/* move the positions of <a>, the analysis of the part of a read of
   length <len> from <from> to <to>, to the whole read */
static void shift_positions(amplicon_t *a, long from, long to, size_t len) {
  if (a->keep_to >= a->keep_from) {
    a->keep_from += from;
    a->keep_to += from;
  }
  for (size_t i = 0; i < a->n_internal; i++) {
    a->internal[i].start += from;
    a->internal[i].end += from;
  }
  /* the primers: in the analyzed orientation (for the minus strand, the
     reverse complement of the read) */
  long shift = a->strand == STRAND_MINUS ? (long)len - 1 - to : from;
  amplicon_end_t *ends[] = {&a->left, &a->right};
  for (size_t e = 0; e < 2; e++)
    if (ends[e]->found) {
      ends[e]->hit.start += shift;
      ends[e]->hit.end += shift;
      ends[e]->hit.terminus_start += shift;
      ends[e]->hit.terminus_end += shift;
    }
}

/* the analysis of a read, with the adapters (see amplicon.h), into
   <out> */
void amplicon_analyze_with_adapters(amplicon_analyzer_t *self,
    const char *seq, size_t len, amplicon_t *out) {
  amplicon_analyze(self, seq, len, self->window, self->window, out);
  out->adapter_start.found = out->adapter_end.found = false;
  if (self->adapters == NULL)
    return;
  /* the adapters, at the ends of the read without a primer */
  bool primer_start = out->strand == STRAND_MINUS ? out->right.found :
                      out->left.found;
  bool primer_end = out->strand == STRAND_MINUS ? out->left.found :
                    out->right.found;
  adapter_hit_t a_start = {.found = false}, a_end = {.found = false};
  if (!primer_start)
    adapters_search_start(self->adapters, seq, len, &a_start);
  if (!primer_end)
    adapters_search_end(self->adapters, seq, len, &a_end);
  if (!a_start.found && !a_end.found)
    return;
  long from = a_start.found ? a_start.end + 1 : 0;
  long to = a_end.found ? a_end.start - 1 : (long)len - 1;
  if (to < from) {
    /* nothing between the adapters */
    amplicon_analyze(self, seq, 0, self->window, self->window, out);
    out->keep_from = 1;
    out->keep_to = 0;
  } else {
    /* the part between the adapters is analyzed again: a primer whose
       outer bases are missing is now at an end, next to the adapter (see
       amplicon.h) */
    amplicon_t whole = *out;
    amplicon_analyze(self, seq + from, (size_t)(to - from + 1),
        a_start.found ? AMPLICON_ADAPTER_WINDOW : self->window,
        a_end.found ? AMPLICON_ADAPTER_WINDOW : self->window, out);
    if (amplicon_n_ends(out) < amplicon_n_ends(&whole))
      *out = whole;
    else
      shift_positions(out, from, to, len);
    if (out->keep_from < from)
      out->keep_from = from;
    if (out->keep_to > to)
      out->keep_to = to;
  }
  out->adapter_start = a_start;
  out->adapter_end = a_end;
}
