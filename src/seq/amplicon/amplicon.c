/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdlib.h>
#include "alloc.h"
#include "amplicon.h"
#include "amplicon-adapters.h"
#include "amplicon-analyzer.h"
#include "amplicon-cut.h"
#include "amplicon-ends.h"
#include "amplicon-locate.h"
#include "fasta.h"
#include "internal-matches.h"
#include "primer-search.h"
#include "reference.h"
#include "scheme.h"
#include "segments.h"
#include "sequtil.h"

amplicon_analyzer_t* amplicon_analyzer_new(const scheme_t *scheme,
    double max_errors, trim_policy_t policy, long window) {
  assert(scheme != NULL);
  assert(window >= 0);
  amplicon_analyzer_t *self;
  XMALLOC(self, sizeof(amplicon_analyzer_t));
  self->scheme = scheme;
  self->max_errors = max_errors;
  self->policy = policy;
  self->window = window;
  self->max_primer_len = 0;
  for (size_t i = 0; i < scheme->n_left; i++)
    if (scheme->left[i].seq_len > self->max_primer_len)
      self->max_primer_len = scheme->left[i].seq_len;
  for (size_t i = 0; i < scheme->n_right; i++)
    if (scheme->right[i].seq_len > self->max_primer_len)
      self->max_primer_len = scheme->right[i].seq_len;
  XMALLOC(self->left_max_errors, sizeof(long) * scheme->n_left);
  XMALLOC(self->right_max_errors, sizeof(long) * scheme->n_right);
  for (size_t i = 0; i < scheme->n_left; i++)
    self->left_max_errors[i] = amplicon_max_errors_of(max_errors,
                                   scheme->left[i].seq_len);
  for (size_t i = 0; i < scheme->n_right; i++)
    self->right_max_errors[i] = amplicon_max_errors_of(max_errors,
                                    scheme->right[i].seq_len);
  size_t stride = self->max_primer_len + 1;
  XMALLOC(self->left_skip, sizeof(long) * stride * scheme->n_left);
  XMALLOC(self->left_edge, sizeof(long) * stride * scheme->n_left);
  XMALLOC(self->right_skip, sizeof(long) * stride * scheme->n_right);
  XMALLOC(self->right_edge, sizeof(long) * stride * scheme->n_right);
  for (size_t i = 0; i < scheme->n_left; i++)
    amplicon_missing_costs(max_errors, scheme->left[i].seq_len,
                           self->left_max_errors[i],
                           self->left_skip + i * stride,
                           self->left_edge + i * stride);
  for (size_t i = 0; i < scheme->n_right; i++)
    amplicon_missing_costs(max_errors, scheme->right[i].seq_len,
                           self->right_max_errors[i],
                           self->right_skip + i * stride,
                           self->right_edge + i * stride);
  self->search = primer_search_new();
  self->internal = internal_matches_new(scheme, self->left_max_errors,
                                        self->right_max_errors);
  self->types = VIRUS_SET_BOTH;
  self->rc = NULL;
  self->rc_cap = 0;
  self->refs = NULL;
  self->ref_segment = NULL;
  self->terminal_window = 0;
  self->probe = NULL;
  self->adapters = NULL;
  self->template = NULL;
  self->template_cap = 0;
  self->skipped = NULL;
  self->skipped_cap = 0;
  return self;
}

void amplicon_analyzer_use_adapters(amplicon_analyzer_t *self,
    adapters_t *adapters) {
  assert(adapters != NULL);
  self->adapters = adapters_any(adapters) ? adapters : NULL;
}

void amplicon_analyzer_use_references(amplicon_analyzer_t *self,
    reference_set_t *refs, long terminal_window) {
  assert(refs != NULL);
  assert(terminal_window >= 0);
  self->refs = refs;
  self->terminal_window = terminal_window;
  XREALLOC(self->probe, reference_set_probe_len(refs) + 1);
  const fasta_t *fa = reference_set_fasta(refs);
  free(self->ref_segment);
  self->ref_segment = NULL;
  if (fa != NULL && fa->n_records > 0) {
    XMALLOC(self->ref_segment, sizeof(int) * fa->n_records);
    for (size_t i = 0; i < fa->n_records; i++) {
      const char *description = fa->records[i].description;
      /* (a reference whose header is a name alone has no segment) */
      self->ref_segment[i] = description == NULL ? 0 :
                             segments_of_description(description);
    }
  }
}

void amplicon_analyzer_restrict_virus(amplicon_analyzer_t *self,
    virus_set_t types) {
  assert(self != NULL);
  self->types = types;
  internal_matches_restrict_virus(self->internal, types);
}

void amplicon_analyzer_restrict_references(amplicon_analyzer_t *self,
    const bool *allowed) {
  if (self->refs != NULL)
    reference_set_restrict(self->refs, allowed);
}

void amplicon_analyzer_delete(amplicon_analyzer_t *self) {
  if (self == NULL)
    return;
  free(self->left_max_errors);
  free(self->right_max_errors);
  free(self->left_skip);
  free(self->left_edge);
  free(self->right_skip);
  free(self->right_edge);
  free(self->ref_segment);
  primer_search_delete(self->search);
  internal_matches_delete(self->internal);
  free(self->rc);
  free(self->probe);
  free(self->template);
  free(self->skipped);
  free(self);
}

/* the number of errors of the primers found (for the comparison of the
   two orientations) */
static long total_errors(const amplicon_t *a) {
  return (a->left.found ? a->left.hit.n_errors : 0) +
         (a->right.found ? a->right.hit.n_errors : 0);
}

/* see amplicon-analyzer.h */
void amplicon_analyze(amplicon_analyzer_t *self, const char *seq, size_t len,
    long start_window, long end_window, amplicon_t *out) {
  amplicon_t plus;
  amplicon_search_layout(self, seq, len, start_window, end_window, &plus);
  amplicon_t minus;
  bool minus_searched = false;
  /*
     The primers of the two ends of a scheme can be similar to each
     other: the tailed primers share their tail, and only their
     conserved termini differ. A primer of one end can therefore be
     found, with errors, at the other end of a read of the opposite
     orientation. Both orientations are thus searched and compared,
     except when the plus strand already gives a complete layout
     without any error, which no other orientation can better.
  */
  if (amplicon_n_ends(&plus) < 2 || total_errors(&plus) > 0) {
    XCHECKSPACEINC(self->rc, &self->rc_cap, len + 1, 1024);
    reverse_complement(seq, len, self->rc);
    self->rc[len] = '\0';
    amplicon_search_layout(self, self->rc, len, end_window, start_window,
                           &minus);
    minus_searched = true;
  }
  /* the orientation with more primers found; with the same number, the
     one with fewer errors; with the same number of errors, the plus
     strand */
  bool use_minus = minus_searched && amplicon_n_ends(&minus) > 0 &&
      (amplicon_n_ends(&minus) > amplicon_n_ends(&plus) ||
       (amplicon_n_ends(&minus) == amplicon_n_ends(&plus) &&
        total_errors(&minus) < total_errors(&plus)));
  *out = use_minus ? minus : plus;
  amplicon_search_internal(self, use_minus ? self->rc : seq, len, out);
  if (self->refs != NULL)
    amplicon_classify_ends(self, use_minus ? self->rc : seq, len, out);
  out->strand = amplicon_n_ends(out) == 0 ? STRAND_NONE :
                (use_minus ? STRAND_MINUS : STRAND_PLUS);
  long from, to;
  amplicon_compute_keep(out, self->policy, len, &from, &to);
  if (minus_searched && to >= from)
    amplicon_cut_facing_primer(out, use_minus ? &plus : &minus, len,
                               &from, &to);
  /* the positions of the original read: for a read of the minus strand,
     the kept part is mirrored back (the read is not reverse-complemented,
     see amplicon.h) */
  if (use_minus) {
    long n = (long)len;
    out->keep_from = to < from ? 1 : n - 1 - to;
    out->keep_to = to < from ? 0 : n - 1 - from;
    for (size_t i = 0; i < out->n_internal; i++) {
      amplicon_match_t *m = out->internal + i;
      long start = n - 1 - m->end;
      m->end = n - 1 - m->start;
      m->start = start;
    }
  } else {
    out->keep_from = from;
    out->keep_to = to;
  }
}

const amplicon_t* amplicon_analyzer_run(amplicon_analyzer_t *self,
    const char *seq, size_t len) {
  assert(self != NULL);
  assert(seq != NULL);
  amplicon_t *out = &self->result;
  amplicon_analyze_with_adapters(self, seq, len, out);
  amplicon_place(self, seq, out);
  amplicon_locate_inner_ends(self, seq, len, out);
  amplicon_find_junction(self, seq, len, out);
  return out;
}
