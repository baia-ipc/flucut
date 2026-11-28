/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "codons.h"
/* the table of the markers, written by the build from
   data/markers.yaml (scripts/markers-table.awk) */
#define MARKERS_TABLE_DEFINE
#include "markers.h"

const marker_def_t* marker_def(size_t marker) {
  assert(marker < MARKER_N_MARKERS);
  return MARKERS_TABLE + marker;
}

/* --- the values read --- */

bool marker_value_of(marker_kind_t kind, const char *nt, size_t window,
    char *out) {
  assert(nt != NULL);
  assert(out != NULL);
  out[0] = '\0';
  char aa[MARKER_WINDOW_NT / 3 + 1];
  if (window < 3 || window > MARKER_WINDOW_NT)
    return false;
  size_t n = codons_translate(nt, window - window % 3, aa);
  switch (kind) {
    case MARKER_RESIDUES:
      if (n > MARKER_VALUE_MAX)
        return false;
      memcpy(out, aa, n + 1);
      return true;
    case MARKER_TO_STOP: {
      const char *stop = strchr(aa, '*');
      /* the last residues before the stop codon, the stop included: the
         end of the protein, read from the stop backwards, so that the
         value is a few codons wide however far into the window the
         protein ends. A window without a stop codon does not say where
         it ends. */
      if (stop == NULL)
        return false;
      size_t before = (size_t)(stop - aa);
      if (before > MARKER_TERMINAL_AA)
        before = MARKER_TERMINAL_AA;
      memcpy(out, stop - before, before);
      out[before] = '*';
      out[before + 1] = '\0';
      return true;
    }
    case MARKER_FRAME_LENGTH: {
      const char *stop = strchr(aa, '*');
      if (stop == NULL)
        return false;
      snprintf(out, MARKER_VALUE_MAX + 1, "%ld", (long)(stop - aa));
      return true;
    }
    default:
      return false;
  }
}

void marker_span_value(long bases, char *out) {
  assert(out != NULL);
  snprintf(out, MARKER_VALUE_MAX + 1, "%ld", bases);
}

bool marker_value_attention(const marker_def_t *def, const char *value,
    const char *reference) {
  assert(def != NULL);
  assert(value != NULL);
  assert(reference != NULL);
  if (def->kind == MARKER_FRAME_LENGTH)
    return def->attention_below > 0 &&
           strtoul(value, NULL, 10) < def->attention_below;
  if (def->kind == MARKER_SPAN) {
    if (def->attention_below == 0 || reference[0] == '\0')
      return false;
    long missing = (long)strtoul(reference, NULL, 10) -
                   (long)strtoul(value, NULL, 10);
    return missing >= (long)def->attention_below;
  }
  for (size_t i = 0; def->attention[i] != '\0' && value[i] != '\0'; i++)
    if (def->attention[i] != '.' && def->attention[i] == value[i])
      return true;
  return false;
}

/* --- the anchors in the references --- */

struct marker_anchors_t {
  size_t n_refs;
  size_t n_found;
  /* one row of MARKER_N_MARKERS per reference */
  marker_reference_t *refs;
};

static marker_reference_t* row(const marker_anchors_t *self, size_t ref) {
  return self->refs + ref * MARKER_N_MARKERS;
}

/*
   The patterns of the table.

   A pattern is a peptide in which a position may be a residue, a '.',
   which any residue matches, or a class "[AS]", which the residues in
   it match. The classes keep an anchor of a position which varies
   between the subtypes specific enough not to be found elsewhere in the
   segment: without them a pattern loose enough to cover every subtype
   also matches sequence which is not the site (docs/markers.md).
*/

/* the residues of <pattern>: a class counts as one */
static size_t pattern_len(const char *pattern) {
  size_t n = 0;
  for (const char *p = pattern; *p != '\0'; p++, n++)
    if (*p == '[')
      while (*p != ']' && *p != '\0')
        p++;
  return n;
}

/* whether the <avail> residues at <aa> begin with <pattern> */
static bool pattern_match(const char *pattern, const char *aa, size_t avail) {
  size_t i = 0;
  for (const char *p = pattern; *p != '\0'; p++, i++) {
    if (i >= avail)
      return false;
    if (*p == '.')
      continue;
    if (*p != '[') {
      if (*p != aa[i])
        return false;
      continue;
    }
    bool in = false;
    while (*++p != ']' && *p != '\0')
      in = in || *p == aa[i];
    if (!in)
      return false;
  }
  return true;
}

/*
   The 0-based position of the first base of the first codon of
   <pattern> in the <len> bases of <seq>, if the pattern occurs exactly
   once in the three reading frames together; -1 if it occurs in none,
   -2 if it occurs more than once. An anchor which could be found in
   more than one place is not used at all: which of them it would find
   in a read is not decided.

   <aa> is a buffer of at least len / 3 + 1 characters.
*/
static long unique_peptide(const char *seq, size_t len, const char *pattern,
    char *aa) {
  long first = -1;
  size_t n = 0;
  size_t plen = pattern_len(pattern);
  for (size_t frame = 0; frame < 3 && len > frame; frame++) {
    size_t got = codons_translate(seq + frame, len - frame - (len - frame) % 3,
                                  aa);
    for (size_t i = 0; i + plen <= got; i++)
      if (pattern_match(pattern, aa + i, got - i)) {
        n++;
        if (first < 0)
          first = (long)(frame + 3 * i);
      }
  }
  return n == 0 ? -1 : (n > 1 ? -2 : first);
}

/*
   Whether the bases of the anchor are found nowhere else in the
   reference within the errors a read is searched with.

   That the peptide occurs once is not enough: the anchor is searched in
   the reads as bases, with errors allowed, so a stretch elsewhere in
   the segment which is near enough to it would be found instead in a
   read whose own site carries an error. The two parts of the reference
   outside the anchor are searched, in both orientations, as a read is.
*/
static bool anchor_alone(const anchor_t *anchor, const char *seq, size_t len,
    long pos, double max_error_rate) {
  anchor_hit_t hit;
  size_t after = (size_t)pos + anchor->len;
  return !anchor_find(anchor, seq, (size_t)pos, max_error_rate, &hit) &&
         !anchor_find(anchor, seq + after, len - after, max_error_rate, &hit);
}

/*
   Build the anchors of <marker> in the reference <r>, whose sequence is
   the <len> bases of <seq>, into <out>; false, and <out> not found, if
   an anchor is not there exactly once or the windows do not fit
*/
static bool anchors_of(marker_reference_t *out, size_t marker, size_t r,
    const char *seq, size_t len, char *aa, double max_error_rate,
    const char **why) {
  const marker_def_t *def = marker_def(marker);
  memset(out, 0, sizeof(marker_reference_t));
  long pos = unique_peptide(seq, len, def->peptide, aa);
  if (pos < 0) {
    *why = pos == -1 ? "the anchor is not in it" :
                       "the anchor is in it more than once";
    return false;
  }
  if (!anchor_at(&out->anchor, r, seq, len, pos, MARKER_ANCHOR_NT)) {
    *why = "the anchor reaches past its end";
    return false;
  }
  if (!anchor_alone(&out->anchor, seq, len, pos, max_error_rate)) {
    *why = "the bases of the anchor are found elsewhere in it as well";
    return false;
  }
  codons_translate(out->anchor.seq, MARKER_ANCHOR_NT, out->anchor_aa);
  if (def->kind == MARKER_SPAN) {
    long pos2 = unique_peptide(seq, len, def->peptide2, aa);
    if (pos2 < 0) {
      *why = pos2 == -1 ? "the second anchor is not in it" :
                          "the second anchor is in it more than once";
      return false;
    }
    if (!anchor_at(&out->anchor2, r, seq, len, pos2, MARKER_ANCHOR_NT)) {
      *why = "the second anchor reaches past its end";
      return false;
    }
    if (!anchor_alone(&out->anchor2, seq, len, pos2, max_error_rate)) {
      *why = "the bases of the second anchor are found elsewhere as well";
      return false;
    }
    marker_span_value(pos2 - pos, out->value);
    out->found = true;
    return true;
  }
  /* the value the reference itself carries */
  long from = def->offset < 0 ? pos + def->offset
                              : pos + (long)MARKER_ANCHOR_NT + def->offset;
  if (from < 0 || (size_t)from + def->window > len) {
    *why = "the window reaches past its end";
    return false;
  }
  if (!marker_value_of(def->kind, seq + from, def->window, out->value)) {
    *why = "the window of the reference gives no value";
    return false;
  }
  out->found = true;
  return true;
}

marker_anchors_t* marker_anchors_new(const fasta_t *references,
    const reference_index_t *index, double max_error_rate, Logger *logger) {
  if (references == NULL || index == NULL || references->n_records == 0)
    return NULL;
  marker_anchors_t *self;
  XCALLOC(self, 1, sizeof(marker_anchors_t));
  self->n_refs = references->n_records;
  XCALLOC(self->refs, self->n_refs * MARKER_N_MARKERS,
          sizeof(marker_reference_t));
  size_t longest = 0;
  for (size_t r = 0; r < self->n_refs; r++)
    if (references->records[r].seq_len > longest)
      longest = references->records[r].seq_len;
  char *aa;
  XMALLOC(aa, longest / 3 + 2);
  for (size_t r = 0; r < self->n_refs; r++) {
    const fasta_record_t *rec = references->records + r;
    int segment = reference_index_segment(index, r);
    marker_reference_t *refs = row(self, r);
    for (size_t m = 0; m < MARKER_N_MARKERS; m++) {
      if (segment < 0 || (size_t)segment != marker_def(m)->segment)
        continue;
      const char *why = NULL;
      if (anchors_of(refs + m, m, r, rec->seq, rec->seq_len, aa,
                     max_error_rate, &why))
        self->n_found++;
      else
        lprintf(logger, LOG_DEBUG, "The marker %s is not read from the "
                "reference %s: %s\n", marker_def(m)->name, rec->name, why);
    }
  }
  free(aa);
  if (self->n_found == 0) {
    marker_anchors_delete(self);
    return NULL;
  }
  return self;
}

void marker_anchors_delete(marker_anchors_t *self) {
  if (self == NULL)
    return;
  free(self->refs);
  free(self);
}

size_t marker_anchors_n(const marker_anchors_t *self) {
  return self == NULL ? 0 : self->n_found;
}

void marker_anchors_chosen(const marker_anchors_t *self, const bool *allowed,
    const marker_reference_t **out) {
  assert(out != NULL);
  for (size_t m = 0; m < MARKER_N_MARKERS; m++)
    out[m] = NULL;
  if (self == NULL)
    return;
  for (size_t r = 0; r < self->n_refs; r++) {
    if (allowed != NULL && !allowed[r])
      continue;
    const marker_reference_t *refs = row(self, r);
    for (size_t m = 0; m < MARKER_N_MARKERS; m++)
      if (out[m] == NULL && refs[m].found)
        out[m] = refs + m;
  }
}
