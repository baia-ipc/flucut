/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "cleavage-site.h"
#include "codons.h"
#include "segments.h"
#include "test.h"

#define HA_SEGMENT 4

/* --- how the residues of a site read --- */

static bool is_basic(char residue) {
  return residue == 'R' || residue == 'K';
}

size_t cleavage_basic_residues(const char *aa) {
  size_t n = 0, len = strlen(aa);
  size_t from = len > 6 ? len - 6 : 0;
  for (size_t i = from; i < len; i++)
    n += is_basic(aa[i]) ? 1 : 0;
  return n;
}

size_t cleavage_reading_of(const char *aa) {
  size_t len = strlen(aa);
  if (len == 0 || strchr(aa, '*') != NULL || strchr(aa, 'X') != NULL)
    return CLEAVAGE_UNREADABLE;
  if (!is_basic(aa[len - 1]))
    return CLEAVAGE_UNREADABLE;
  return cleavage_basic_residues(aa) >= 4 ? CLEAVAGE_MULTIBASIC
                                          : CLEAVAGE_MONOBASIC;
}

const char* cleavage_reading_name(size_t reading) {
  static const char *const names[CLEAVAGE_N_READINGS] = {
    "multibasic", "monobasic", "unreadable"};
  assert(reading < CLEAVAGE_N_READINGS);
  return names[reading];
}

struct cleavage_anchors_t {
  size_t n_refs;
  size_t n_found;
  /* one per reference; the anchor is found only for the HA ones */
  cleavage_reference_t *refs;
};

cleavage_anchors_t* cleavage_anchors_new(const fasta_t *references,
    const reference_index_t *index, Logger *logger) {
  if (references == NULL || index == NULL || references->n_records == 0)
    return NULL;
  cleavage_anchors_t *self;
  XCALLOC(self, 1, sizeof(cleavage_anchors_t));
  self->n_refs = references->n_records;
  XCALLOC(self->refs, self->n_refs, sizeof(cleavage_reference_t));
  for (size_t r = 0; r < self->n_refs; r++) {
    if (reference_index_segment(index, r) != HA_SEGMENT)
      continue;
    const fasta_record_t *rec = references->records + r;
    cleavage_reference_t *cr = self->refs + r;
    /* what the anchor is worth in this reference, reported and not
       acted on: how often its peptide occurs, and whether its bases are
       found elsewhere in the segment within the errors a read is
       searched with. A second occurrence of either would mean the
       anchor could be found in the wrong place in a read */
    long first = -1;
    size_t n_peptide = anchor_peptide_occurrences(rec->seq, rec->seq_len,
        CLEAVAGE_FUSION_PEPTIDE, &first);
    if (n_peptide != 1)
      lprintf(logger, LOG_DEBUG, "The fusion peptide occurs %zu times in "
              "the HA reference %s\n", n_peptide, rec->name);
    if (!anchor_of_peptide(&cr->anchor, r, rec->seq, rec->seq_len,
                           CLEAVAGE_FUSION_PEPTIDE, CLEAVAGE_ANCHOR_NT)) {
      lprintf(logger, LOG_DEBUG,
              "No fusion peptide in the HA reference %s: the cleavage "
              "site is not read from it\n", rec->name);
      continue;
    }
    /* the site the reference itself carries, which the reads are
       compared with; a reference whose HA2 begins too near its start
       has none */
    if (cr->anchor.ref_pos < CLEAVAGE_SITE_NT) {
      cr->anchor.found = false;
      continue;
    }
    memcpy(cr->nt, rec->seq + cr->anchor.ref_pos - CLEAVAGE_SITE_NT,
           CLEAVAGE_SITE_NT);
    cr->nt[CLEAVAGE_SITE_NT] = '\0';
    codons_translate(cr->nt, CLEAVAGE_SITE_NT, cr->aa);
    cr->reading = cleavage_reading_of(cr->aa);
    codons_translate(cr->anchor.seq, CLEAVAGE_ANCHOR_NT, cr->fusion);
    if (!anchor_bases_alone(&cr->anchor, rec->seq, rec->seq_len, first,
                            CLEAVAGE_ANCHOR_ERROR_RATE))
      lprintf(logger, LOG_DEBUG, "The bases of the fusion peptide anchor "
              "are found elsewhere in the HA reference %s as well\n",
              rec->name);
    self->n_found++;
  }
  if (self->n_found == 0) {
    cleavage_anchors_delete(self);
    return NULL;
  }
  return self;
}

const cleavage_reference_t* cleavage_anchors_chosen(
    const cleavage_anchors_t *self, const bool *allowed) {
  if (self == NULL)
    return NULL;
  for (size_t r = 0; r < self->n_refs; r++)
    if (self->refs[r].anchor.found && (allowed == NULL || allowed[r]))
      return self->refs + r;
  return NULL;
}

size_t cleavage_anchors_n(const cleavage_anchors_t *self) {
  return self == NULL ? 0 : self->n_found;
}

void cleavage_anchors_delete(cleavage_anchors_t *self) {
  if (self == NULL)
    return;
  free(self->refs);
  free(self);
}
