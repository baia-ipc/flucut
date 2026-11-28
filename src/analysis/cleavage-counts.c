/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "cleavage-counts.h"
#include "codons.h"
#include "test.h"

/* --- the motifs counted --- */

void cleavage_counts_clear(cleavage_counts_t *self) {
  assert(self != NULL);
  memset(self, 0, sizeof(cleavage_counts_t));
}

void cleavage_counts_set_reference(cleavage_counts_t *self, const char *aa,
    const char *fusion) {
  assert(self != NULL);
  assert(aa != NULL);
  assert(fusion != NULL);
  if (self->reference_aa[0] != '\0')
    return;
  snprintf(self->reference_aa, sizeof(self->reference_aa), "%s", aa);
  snprintf(self->reference_fusion, sizeof(self->reference_fusion), "%s",
           fusion);
}

void cleavage_counts_searched(cleavage_counts_t *self) {
  assert(self != NULL);
  self->searched++;
}

/* FNV-1a over the bases of a motif */
static size_t hash(const char *nt) {
  size_t h = 2166136261u;
  for (size_t i = 0; i < CLEAVAGE_SITE_NT; i++) {
    h ^= (unsigned char)nt[i];
    h *= 16777619u;
  }
  return h;
}

/* the slot of <nt> in the table: the one holding it, or the first empty
   one after its place */
static size_t slot_of(const cleavage_counts_t *self, const char *nt) {
  size_t s = hash(nt) & (CLEAVAGE_SLOTS - 1);
  while (self->slots[s] != 0 &&
         memcmp(self->motifs[self->slots[s] - 1].nt, nt,
                CLEAVAGE_SITE_NT) != 0)
    s = (s + 1) & (CLEAVAGE_SLOTS - 1);
  return s;
}

/* add <reads> reads of the motif <nt> */
static void add_reads(cleavage_counts_t *self, const char *nt,
    size_t reads) {
  size_t s = slot_of(self, nt);
  if (self->slots[s] != 0) {
    self->motifs[self->slots[s] - 1].reads += reads;
    return;
  }
  if (self->n_motifs == CLEAVAGE_MAX_MOTIFS) {
    self->overflow += reads;
    return;
  }
  memcpy(self->motifs[self->n_motifs].nt, nt, CLEAVAGE_SITE_NT);
  self->motifs[self->n_motifs].nt[CLEAVAGE_SITE_NT] = '\0';
  self->motifs[self->n_motifs].reads = reads;
  self->n_motifs++;
  self->slots[s] = self->n_motifs;
}

void cleavage_counts_add(cleavage_counts_t *self, const char *nt) {
  assert(self != NULL);
  assert(nt != NULL);
  add_reads(self, nt, 1);
}

/* the slot of <aa> in the table of a control window */
static size_t control_slot(const cleavage_tally_t *c, const char *aa) {
  size_t h = 2166136261u;
  for (size_t i = 0; i < CLEAVAGE_SITE_AA; i++) {
    h ^= (unsigned char)aa[i];
    h *= 16777619u;
  }
  size_t s = h & (CLEAVAGE_CONTROL_SLOTS - 1);
  while (c->slots[s] != 0 &&
         memcmp(c->entries[c->slots[s] - 1].aa, aa, CLEAVAGE_SITE_AA) != 0)
    s = (s + 1) & (CLEAVAGE_CONTROL_SLOTS - 1);
  return s;
}

static void control_add(cleavage_tally_t *c, const char *aa,
    size_t reads) {
  c->reads += reads;
  size_t s = control_slot(c, aa);
  if (c->slots[s] != 0) {
    c->entries[c->slots[s] - 1].reads += reads;
    return;
  }
  if (c->n == CLEAVAGE_CONTROL_MOTIFS) {
    c->overflow += reads;
    return;
  }
  memcpy(c->entries[c->n].aa, aa, CLEAVAGE_SITE_AA);
  c->entries[c->n].aa[CLEAVAGE_SITE_AA] = '\0';
  c->entries[c->n].reads = reads;
  c->n++;
  c->slots[s] = c->n;
}

void cleavage_counts_add_control(cleavage_counts_t *self, size_t i,
    const char *nt) {
  assert(self != NULL);
  assert(nt != NULL);
  assert(i < CLEAVAGE_N_CONTROLS);
  char aa[CLEAVAGE_SITE_AA + 1];
  codons_translate(nt, CLEAVAGE_SITE_NT, aa);
  control_add(self->control + i, aa, 1);
}

void cleavage_counts_add_fusion(cleavage_counts_t *self, const char *nt) {
  assert(self != NULL);
  assert(nt != NULL);
  char aa[CLEAVAGE_SITE_AA + 1];
  codons_translate(nt, CLEAVAGE_ANCHOR_NT, aa);
  control_add(&self->fusion, aa, 1);
}

size_t cleavage_counts_reads(const cleavage_counts_t *self) {
  assert(self != NULL);
  size_t n = self->overflow;
  for (size_t i = 0; i < self->n_motifs; i++)
    n += self->motifs[i].reads;
  return n;
}

void cleavage_counts_merge(cleavage_counts_t *self,
    const cleavage_counts_t *other) {
  assert(self != NULL);
  assert(other != NULL);
  self->searched += other->searched;
  self->overflow += other->overflow;
  for (size_t i = 0; i < CLEAVAGE_N_CONTROLS; i++) {
    self->control[i].overflow += other->control[i].overflow;
    self->control[i].reads += other->control[i].overflow;
    for (size_t k = 0; k < other->control[i].n; k++)
      control_add(self->control + i, other->control[i].entries[k].aa,
                  other->control[i].entries[k].reads);
  }
  if (other->reference_aa[0] != '\0')
    cleavage_counts_set_reference(self, other->reference_aa,
                                  other->reference_fusion);
  self->fusion.overflow += other->fusion.overflow;
  self->fusion.reads += other->fusion.overflow;
  for (size_t k = 0; k < other->fusion.n; k++)
    control_add(&self->fusion, other->fusion.entries[k].aa,
                other->fusion.entries[k].reads);
  for (size_t i = 0; i < other->n_motifs; i++)
    add_reads(self, other->motifs[i].nt, other->motifs[i].reads);
}

void cleavage_counts_save(const cleavage_counts_t *self, FILE *file) {
  assert(self != NULL);
  assert(file != NULL);
  fprintf(file, "cleavage-reads %zu %zu\n", self->searched, self->overflow);
  for (size_t i = 0; i < CLEAVAGE_N_CONTROLS; i++) {
    if (self->control[i].overflow > 0)
      fprintf(file, "cleavage-control-over %zu %zu\n", i,
              self->control[i].overflow);
    for (size_t k = 0; k < self->control[i].n; k++)
      fprintf(file, "cleavage-control %zu %s %zu\n", i,
              self->control[i].entries[k].aa,
              self->control[i].entries[k].reads);
  }
  if (self->reference_aa[0] != '\0')
    fprintf(file, "cleavage-reference %s %s\n", self->reference_aa,
            self->reference_fusion);
  if (self->fusion.overflow > 0)
    fprintf(file, "cleavage-fusion-over %zu\n", self->fusion.overflow);
  for (size_t k = 0; k < self->fusion.n; k++)
    fprintf(file, "cleavage-fusion %s %zu\n", self->fusion.entries[k].aa,
            self->fusion.entries[k].reads);
  for (size_t i = 0; i < self->n_motifs; i++)
    fprintf(file, "cleavage %s %zu\n", self->motifs[i].nt,
            self->motifs[i].reads);
}

bool cleavage_counts_load_line(cleavage_counts_t *self, char *line) {
  assert(self != NULL);
  assert(line != NULL);
  if (strncmp(line, "cleavage-reads ", 15) == 0) {
    char *rest;
    self->searched += strtoull(line + 15, &rest, 10);
    self->overflow += strtoull(rest, NULL, 10);
    return true;
  }
  if (strncmp(line, "cleavage-control-over ", 22) == 0) {
    char *rest;
    size_t i = strtoull(line + 22, &rest, 10);
    if (i >= CLEAVAGE_N_CONTROLS)
      return false;
    size_t n = strtoull(rest, NULL, 10);
    self->control[i].overflow += n;
    self->control[i].reads += n;
    return true;
  }
  if (strncmp(line, "cleavage-control ", 17) == 0) {
    char *rest;
    size_t i = strtoull(line + 17, &rest, 10);
    char aa[CLEAVAGE_SITE_AA + 2];
    if (i >= CLEAVAGE_N_CONTROLS || sscanf(rest, "%9s", aa) != 1 ||
        strlen(aa) != CLEAVAGE_SITE_AA)
      return false;
    rest = strchr(rest + 1, ' ');
    if (rest == NULL)
      return false;
    control_add(self->control + i, aa, strtoull(rest, NULL, 10));
    return true;
  }
  if (strncmp(line, "cleavage-fusion-over ", 21) == 0) {
    size_t n = strtoull(line + 21, NULL, 10);
    self->fusion.overflow += n;
    self->fusion.reads += n;
    return true;
  }
  if (strncmp(line, "cleavage-fusion ", 16) == 0) {
    char aa[CLEAVAGE_SITE_AA + 2];
    char *rest;
    if (sscanf(line + 16, "%9s", aa) != 1 ||
        strlen(aa) != CLEAVAGE_SITE_AA ||
        (rest = strchr(line + 16, ' ')) == NULL)
      return false;
    control_add(&self->fusion, aa, strtoull(rest, NULL, 10));
    return true;
  }
  if (strncmp(line, "cleavage-reference ", 19) == 0) {
    char aa[CLEAVAGE_SITE_AA + 2], fusion[CLEAVAGE_SITE_AA + 2] = "";
    if (sscanf(line + 19, "%9s %9s", aa, fusion) < 1)
      return false;
    cleavage_counts_set_reference(self, aa, fusion);
    return true;
  }
  if (strncmp(line, "cleavage ", 9) != 0)
    return false;
  char *nt = line + 9;
  char *sep = strchr(nt, ' ');
  if (sep == NULL || sep - nt != CLEAVAGE_SITE_NT)
    return false;
  char motif[CLEAVAGE_SITE_NT + 1];
  memcpy(motif, nt, CLEAVAGE_SITE_NT);
  motif[CLEAVAGE_SITE_NT] = '\0';
  add_reads(self, motif, strtoull(sep + 1, NULL, 10));
  return true;
}

/* --- what the reads say --- */
/* the readings carried by enough of the reads to be more than errors */
static size_t supported_readings(const cleavage_motif_t *proteins, size_t n,
    size_t reads) {
  size_t supported = 0;
  for (size_t i = 0; i < n; i++)
    if (proteins[i].reads >= CLEAVAGE_MIN_SITE_READS &&
        (double)proteins[i].reads >= CLEAVAGE_MIN_SITE_SHARE * (double)reads)
      supported++;
  return supported;
}

void cleavage_result_of(const cleavage_counts_t *counts, size_t chemistry,
    cleavage_result_t *out) {
  assert(counts != NULL);
  assert(out != NULL);
  memset(out, 0, sizeof(cleavage_result_t));
  out->chemistry = chemistry;
  out->reference_reading = CLEAVAGE_N_READINGS;
  snprintf(out->reference_aa, sizeof(out->reference_aa), "%s",
           counts->reference_aa);
  if (out->reference_aa[0] != '\0')
    out->reference_reading = cleavage_reading_of(out->reference_aa);
  out->searched = counts->searched;
  out->reads = cleavage_counts_reads(counts);
  out->distinct = counts->n_motifs;
  out->reading = CLEAVAGE_UNREADABLE;
  out->verdict = CLEAVAGE_TOO_FEW_READS;
  if (out->reads == 0 || counts->n_motifs == 0)
    return;
  /* The reads are gathered by the residues they carry, not by the bases:
     a read with an error in the third base of a codon usually still
     reads as the same site, and it is the residues which decide the
     pathotype. The spellings are few enough (CLEAVAGE_MAX_MOTIFS) to
     gather them exactly; this runs once for a report, not for a read. */
  cleavage_motif_t *proteins;
  XCALLOC(proteins, counts->n_motifs, sizeof(cleavage_motif_t));
  size_t *spelling_reads;
  XCALLOC(spelling_reads, counts->n_motifs, sizeof(size_t));
  size_t n = 0, multibasic = 0;
  char aa[CLEAVAGE_SITE_AA + 1];
  for (size_t i = 0; i < counts->n_motifs; i++) {
    codons_translate(counts->motifs[i].nt, CLEAVAGE_SITE_NT, aa);
    size_t reading = cleavage_reading_of(aa);
    if (reading == CLEAVAGE_MULTIBASIC)
      multibasic += counts->motifs[i].reads;
    size_t k = 0;
    while (k < n && strcmp(proteins[k].aa, aa) != 0)
      k++;
    if (k == n) {
      memcpy(proteins[k].aa, aa, CLEAVAGE_SITE_AA + 1);
      proteins[k].reading = reading;
      n++;
    }
    proteins[k].reads += counts->motifs[i].reads;
    /* the most frequent spelling of these residues; among equally
       frequent ones the first in alphabetical order */
    if (counts->motifs[i].reads > spelling_reads[k] ||
        (counts->motifs[i].reads == spelling_reads[k] &&
         strcmp(counts->motifs[i].nt, proteins[k].nt) < 0)) {
      memcpy(proteins[k].nt, counts->motifs[i].nt, CLEAVAGE_SITE_NT + 1);
      spelling_reads[k] = counts->motifs[i].reads;
    }
  }
  out->distinct_proteins = n;
  out->supported_proteins = supported_readings(proteins, n, out->reads);
  /* the most frequent ones, most frequent first; among equally frequent
     residues the first in alphabetical order, so that the result does
     not depend on the order of the reads */
  out->n_top = n < CLEAVAGE_TOP_MOTIFS ? n : CLEAVAGE_TOP_MOTIFS;
  for (size_t i = 0; i < out->n_top; i++) {
    size_t best = i;
    for (size_t j = i + 1; j < n; j++)
      if (proteins[j].reads > proteins[best].reads ||
          (proteins[j].reads == proteins[best].reads &&
           strcmp(proteins[j].aa, proteins[best].aa) < 0))
        best = j;
    cleavage_motif_t t = proteins[i];
    proteins[i] = proteins[best];
    proteins[best] = t;
    out->top[i] = proteins[i];
  }
  free(proteins);
  free(spelling_reads);
  out->found = true;
  memcpy(out->aa, out->top[0].aa, CLEAVAGE_SITE_AA + 1);
  memcpy(out->nt, out->top[0].nt, CLEAVAGE_SITE_NT + 1);
  out->motif_reads = out->top[0].reads;
  out->support = (double)out->motif_reads / (double)out->reads;
  out->basic = cleavage_basic_residues(out->aa);
  out->reading = out->top[0].reading;
  /* the share of the reads whose own site is multibasic, whatever its
     sequence: a lower bound, as the reads over the overflow are not
     counted */
  out->multibasic_share = (double)multibasic / (double)out->reads;
  /* the control windows, pooled: each voted on separately, the reads
     agreeing with its own majority added up */
  size_t agree = 0, total = 0;
  for (size_t i = 0; i < CLEAVAGE_N_CONTROLS; i++) {
    const cleavage_tally_t *c = counts->control + i;
    size_t best = 0;
    for (size_t k = 1; k < c->n; k++)
      if (c->entries[k].reads > c->entries[best].reads)
        best = k;
    if (c->n > 0)
      agree += c->entries[best].reads;
    total += c->reads;
  }
  out->control_reads = total;
  out->control_support = total > 0 ? (double)agree / (double)total : 0.0;
  snprintf(out->reference_fusion, sizeof(out->reference_fusion), "%s",
           counts->reference_fusion);
  const cleavage_tally_t *f = &counts->fusion;
  if (f->reads > 0 && f->n > 0) {
    size_t best = 0;
    for (size_t k = 1; k < f->n; k++)
      if (f->entries[k].reads > f->entries[best].reads ||
          (f->entries[k].reads == f->entries[best].reads &&
           strcmp(f->entries[k].aa, f->entries[best].aa) < 0))
        best = k;
    snprintf(out->fusion_aa, sizeof(out->fusion_aa), "%s",
             f->entries[best].aa);
    out->fusion_reads = f->reads;
    out->fusion_support = (double)f->entries[best].reads / (double)f->reads;
  }
  /* the reads of a sample which holds more than one haemagglutinin do
     not read like its own reference; the motif reported is then not
     necessarily the one of its own virus */
  out->disagrees_with_reference =
      out->reference_reading < CLEAVAGE_N_READINGS &&
      out->reference_reading != out->reading;
  out->control_weak = out->control_reads > 0 &&
      out->control_support < CLEAVAGE_MIN_CONTROL_SUPPORT;
  if (out->reads < CLEAVAGE_MIN_READS)
    out->verdict = CLEAVAGE_TOO_FEW_READS;
  else if (out->disagrees_with_reference || out->control_weak ||
           out->support < CLEAVAGE_MIN_SUPPORT)
    out->verdict = CLEAVAGE_CONFIRM;
  else
    out->verdict = CLEAVAGE_FIRM;
}

/* whether <model> names the chemistry <tag> (e.g. "r10"): the basecalling
   models are named dna_r10.4.1_e8.2_400bps_hac@v4.2.0 and the like */
static bool model_has(const char *model, const char *tag) {
  const char *p = strstr(model, tag);
  /* the tag must be a word of the name, not part of a longer one */
  return p != NULL && (p == model || !((p[-1] >= '0' && p[-1] <= '9') ||
                                       (p[-1] >= 'a' && p[-1] <= 'z')));
}

size_t cleavage_chemistry(bool illumina, const char *basecall_model) {
  if (illumina)
    return CLEAVAGE_CHEMISTRY_ILLUMINA;
  if (basecall_model == NULL || basecall_model[0] == '\0')
    return CLEAVAGE_CHEMISTRY_UNKNOWN;
  if (model_has(basecall_model, "r10"))
    return CLEAVAGE_CHEMISTRY_R10;
  if (model_has(basecall_model, "r9"))
    return CLEAVAGE_CHEMISTRY_R9;
  return CLEAVAGE_CHEMISTRY_UNKNOWN;
}

const char* cleavage_chemistry_name(size_t chemistry) {
  static const char *const names[CLEAVAGE_N_CHEMISTRIES] = {
    "unknown", "R9", "R10", "Illumina"};
  assert(chemistry < CLEAVAGE_N_CHEMISTRIES);
  return names[chemistry];
}

const char* cleavage_verdict_name(size_t verdict) {
  static const char *const names[CLEAVAGE_N_VERDICTS] = {
    "firm", "confirm", "too few reads"};
  assert(verdict < CLEAVAGE_N_VERDICTS);
  return names[verdict];
}
