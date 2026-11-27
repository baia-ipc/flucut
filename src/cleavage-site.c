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

static size_t reading_of(const char *aa);

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
    cr->reading = reading_of(cr->aa);
    codons_translate(cr->anchor.seq, CLEAVAGE_ANCHOR_NT, cr->fusion);
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
    for (size_t k = 0; k < other->control[i].n; k++)
      control_add(self->control + i, other->control[i].entries[k].aa,
                  other->control[i].entries[k].reads);
  }
  if (other->reference_aa[0] != '\0')
    cleavage_counts_set_reference(self, other->reference_aa,
                                  other->reference_fusion);
  self->fusion.overflow += other->fusion.overflow;
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
    self->control[i].overflow += strtoull(rest, NULL, 10);
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
    self->fusion.overflow += strtoull(line + 21, NULL, 10);
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

static bool is_basic(char residue) {
  return residue == 'R' || residue == 'K';
}

/* the basic residues among the last 6 of the HA1 end */
static size_t basic_run(const char *aa) {
  size_t n = 0, len = strlen(aa);
  size_t from = len > 6 ? len - 6 : 0;
  for (size_t i = from; i < len; i++)
    n += is_basic(aa[i]) ? 1 : 0;
  return n;
}

static size_t reading_of(const char *aa) {
  size_t len = strlen(aa);
  if (len == 0 || strchr(aa, '*') != NULL || strchr(aa, 'X') != NULL)
    return CLEAVAGE_UNREADABLE;
  if (!is_basic(aa[len - 1]))
    return CLEAVAGE_UNREADABLE;
  return basic_run(aa) >= 4 ? CLEAVAGE_MULTIBASIC : CLEAVAGE_MONOBASIC;
}

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
    out->reference_reading = reading_of(out->reference_aa);
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
    size_t reading = reading_of(aa);
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
    /* the most frequent spelling of these residues */
    if (counts->motifs[i].reads > spelling_reads[k]) {
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
  out->basic = basic_run(out->aa);
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
      if (f->entries[k].reads > f->entries[best].reads)
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

const char* cleavage_reading_name(size_t reading) {
  static const char *const names[CLEAVAGE_N_READINGS] = {
    "multibasic", "monobasic", "unreadable"};
  assert(reading < CLEAVAGE_N_READINGS);
  return names[reading];
}

const char* cleavage_verdict_name(size_t verdict) {
  static const char *const names[CLEAVAGE_N_VERDICTS] = {
    "firm", "confirm", "too few reads"};
  assert(verdict < CLEAVAGE_N_VERDICTS);
  return names[verdict];
}

/* --- unit tests --- */

/* the cleavage site of an H5 of clade 2.3.4.4b (PQRERRRKR) and of a low
   pathogenic H9 (PARSSR), as the 8 codons before HA2 */
#define HPAI_NT "CAAAGAGAGAGAAGAAGAAAAAGA"
#define LPAI_NT "AATGTTCCTGCAAGATCTTCAAGA"

static void add_n(cleavage_counts_t *c, const char *nt, size_t n) {
  for (size_t i = 0; i < n; i++) {
    cleavage_counts_searched(c);
    cleavage_counts_add(c, nt);
  }
}

static void test_counts(void) {
  cleavage_counts_t c;
  cleavage_counts_clear(&c);
  test_eq_size(cleavage_counts_reads(&c), 0);
  cleavage_result_t r;
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_assert(!r.found);
  test_eq_size(r.verdict, CLEAVAGE_TOO_FEW_READS);

  /* 70 reads of the motif, 20 of a spelling of the same residues (the
     last codon AGG instead of AGA, both R) and 10 whose last residue is
     an isoleucine instead */
  add_n(&c, HPAI_NT, 70);
  add_n(&c, "CAAAGAGAGAGAAGAAGAAAAAGG", 20);
  add_n(&c, "CAAAGAGAGAGAAGAAGAAAAATA", 10);
  test_eq_size(cleavage_counts_reads(&c), 100);
  test_eq_size(c.n_motifs, 3);
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_assert(r.found);
  test_eq_str(r.aa, "QRERRRKR");
  /* the most frequent spelling of those residues */
  test_eq_str(r.nt, HPAI_NT);
  test_eq_size(r.reads, 100);
  /* the two spellings of the same residues are counted together */
  test_eq_size(r.motif_reads, 90);
  test_eq_size(r.distinct, 3);
  test_eq_size(r.distinct_proteins, 2);
  /* both readings are carried by at least a tenth of the reads */
  test_eq_size(r.supported_proteins, 2);
  test_eq_size(r.basic, 5);
  test_eq_size(r.reading, CLEAVAGE_MULTIBASIC);
  test_eq_size(r.verdict, CLEAVAGE_FIRM);
  test_assert(r.support > 0.89 && r.support < 0.91);
  test_assert(r.multibasic_share > 0.89 && r.multibasic_share < 0.91);
  /* without a reference site nothing is compared with it */
  test_eq_size(r.reference_reading, CLEAVAGE_N_READINGS);
  test_assert(!r.disagrees_with_reference);
  /* the sites found, most frequent first */
  test_eq_size(r.n_top, 2);
  test_eq_str(r.top[0].aa, "QRERRRKR");
  test_eq_size(r.top[0].reads, 90);
  test_eq_size(r.top[0].reading, CLEAVAGE_MULTIBASIC);
  test_eq_str(r.top[1].aa, "QRERRRKI");
  test_eq_size(r.top[1].reads, 10);
  test_eq_size(r.top[1].reading, CLEAVAGE_UNREADABLE);

  /* the reads read like the reference: still firm */
  cleavage_counts_set_reference(&c, "LREKRRKR", "GLFGAIAG");
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_str(r.reference_aa, "LREKRRKR");
  test_eq_size(r.reference_reading, CLEAVAGE_MULTIBASIC);
  test_assert(!r.disagrees_with_reference);
  test_eq_size(r.verdict, CLEAVAGE_FIRM);
  /* a sample whose reference is multibasic but whose reads read as a
     low pathogenic site holds the reads of another haemagglutinin: it
     is never firm, however many reads agree */
  cleavage_counts_t mixed;
  cleavage_counts_clear(&mixed);
  cleavage_counts_set_reference(&mixed, "LREKRRKR", "GLFGAIAG");
  add_n(&mixed, LPAI_NT, 100);
  cleavage_result_of(&mixed, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_str(r.aa, "NVPARSSR");
  test_assert(r.support > 0.99);
  test_assert(r.disagrees_with_reference);
  test_eq_size(r.verdict, CLEAVAGE_CONFIRM);

  /* the verdict does not depend on the chemistry: the same reads are
     judged the same whatever they were sequenced with */
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R9, &r);
  test_eq_size(r.verdict, CLEAVAGE_FIRM);
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_ILLUMINA, &r);
  test_eq_size(r.verdict, CLEAVAGE_FIRM);

  /* too few reads */
  cleavage_counts_t few;
  cleavage_counts_clear(&few);
  add_n(&few, HPAI_NT, CLEAVAGE_MIN_READS - 1);
  cleavage_result_of(&few, CLEAVAGE_CHEMISTRY_R10, &r);
  test_assert(r.found);
  test_eq_size(r.verdict, CLEAVAGE_TOO_FEW_READS);
  add_n(&few, HPAI_NT, 1);
  cleavage_result_of(&few, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_size(r.verdict, CLEAVAGE_FIRM);

  /* merging adds the counts of both */
  cleavage_counts_t other;
  cleavage_counts_clear(&other);
  add_n(&other, HPAI_NT, 5);
  add_n(&other, LPAI_NT, 1);
  cleavage_counts_merge(&c, &other);
  test_eq_size(cleavage_counts_reads(&c), 106);
  test_eq_size(c.searched, 106);
  test_eq_size(c.n_motifs, 4);
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_size(r.motif_reads, 95);
}

/* a reading of the site carried by many reads and by a share of them,
   and readings carried by one read each: the raw counts include the
   single reads, the count of supported readings does not */
static void test_supported_readings(void) {
  cleavage_counts_t c;
  cleavage_counts_clear(&c);
  add_n(&c, HPAI_NT, 60);
  add_n(&c, LPAI_NT, 30);
  /* ten reads with a base of their own changed in the first 10 codon
     positions after the first: ten spellings, one read each */
  char err[CLEAVAGE_SITE_NT + 1];
  for (size_t i = 0; i < 10; i++) {
    memcpy(err, HPAI_NT, sizeof(err));
    err[3 + i] = err[3 + i] == 'T' ? 'C' : 'T';
    add_n(&c, err, 1);
  }
  cleavage_result_t r;
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_size(r.reads, 100);
  test_eq_size(r.distinct, 12);
  /* the readings of the single reads are counted here ... */
  test_assert(r.distinct_proteins > 3);
  /* ... and not here: the two readings of the sample */
  test_eq_size(r.supported_proteins, 2);
  test_eq_str(r.top[0].aa, "QRERRRKR");
  test_eq_size(r.top[1].reads, 30);

  /* a reading carried by two reads, below the share, is not supported;
     at the share it is */
  cleavage_counts_t d;
  cleavage_counts_clear(&d);
  add_n(&d, HPAI_NT, 98);
  add_n(&d, LPAI_NT, 2);
  cleavage_result_of(&d, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_size(r.supported_proteins, 1);
  add_n(&d, LPAI_NT, 10);
  cleavage_result_of(&d, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_size(r.supported_proteins, 2);
}

/* the control windows: reads which agree at the site but nowhere else
   are not firm */
static void test_control(void) {
  cleavage_counts_t c;
  cleavage_counts_clear(&c);
  add_n(&c, HPAI_NT, 100);
  cleavage_result_t r;
  /* no control window read: nothing is known about the rest, and the
     site alone decides */
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_size(r.control_reads, 0);
  test_assert(!r.control_weak);
  test_eq_size(r.verdict, CLEAVAGE_FIRM);
  /* control windows the reads agree on (the second motif differs in a
     residue, not only in a base: the control counts residues too) */
  for (size_t i = 0; i < 90; i++)
    for (size_t w = 0; w < CLEAVAGE_N_CONTROLS; w++)
      cleavage_counts_add_control(&c, w, "TTTATAGAGGGAGGATGGCAGGGA");
  for (size_t i = 0; i < 10; i++)
    for (size_t w = 0; w < CLEAVAGE_N_CONTROLS; w++)
      cleavage_counts_add_control(&c, w, "TTTATAGAGGGAGGATGGCAGTCA");
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_size(r.control_reads, 200);
  test_assert(r.control_support > 0.89 && r.control_support < 0.91);
  test_assert(!r.control_weak);
  test_eq_size(r.verdict, CLEAVAGE_FIRM);
  /* reads which agree nowhere, not even where nothing is difficult:
     the site is not wrong so much as unknown */
  cleavage_counts_t noisy;
  cleavage_counts_clear(&noisy);
  add_n(&noisy, HPAI_NT, 100);
  char nt[CLEAVAGE_SITE_NT + 1];
  memcpy(nt, "TTTATAGAGGGAGGATGGCAGGGA", sizeof(nt));
  for (size_t i = 0; i < 100; i++) {
    for (size_t j = 0; j < 4; j++)
      nt[j] = "ACGT"[(i >> (2 * j)) & 3];
    for (size_t w = 0; w < CLEAVAGE_N_CONTROLS; w++)
      cleavage_counts_add_control(&noisy, w, nt);
  }
  cleavage_result_of(&noisy, CLEAVAGE_CHEMISTRY_R10, &r);
  test_assert(r.control_support < CLEAVAGE_MIN_CONTROL_SUPPORT);
  test_assert(r.control_weak);
  test_eq_size(r.verdict, CLEAVAGE_CONFIRM);
  /* the site itself is still reported, with its own support */
  test_eq_str(r.aa, "QRERRRKR");
  test_assert(r.support > 0.99);
}

/* the fusion peptide: the anchor as the reads carry it */
static void test_fusion(void) {
  cleavage_counts_t c;
  cleavage_counts_clear(&c);
  cleavage_counts_set_reference(&c, "LREKRRKR", "GLFGAIAG");
  add_n(&c, HPAI_NT, 100);
  cleavage_result_t r;
  /* no read read at the anchor yet */
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_str(r.reference_fusion, "GLFGAIAG");
  test_eq_str(r.fusion_aa, "");
  test_eq_size(r.fusion_reads, 0);
  /* 80 reads carry the peptide of the reference, 20 another */
  for (size_t i = 0; i < 80; i++)
    cleavage_counts_add_fusion(&c, "GGACTGTTTGGGGCTATAGCAGGC");
  for (size_t i = 0; i < 20; i++)
    cleavage_counts_add_fusion(&c, "GGAATTTTTGGGGCTATAGCAGGC");
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_str(r.fusion_aa, "GLFGAIAG");
  test_eq_size(r.fusion_reads, 100);
  test_assert(r.fusion_support > 0.79 && r.fusion_support < 0.81);
  /* the other peptide is the one of the H3 of horses and pigs */
  cleavage_counts_t h3;
  cleavage_counts_clear(&h3);
  for (size_t i = 0; i < 50; i++)
    cleavage_counts_add_fusion(&h3, "GGAATTTTTGGGGCTATAGCAGGC");
  add_n(&h3, HPAI_NT, 50);
  cleavage_result_of(&h3, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_str(r.fusion_aa, "GIFGAIAG");
  test_assert(r.fusion_support > 0.99);
}

static void test_readings(void) {
  cleavage_counts_t c;
  cleavage_result_t r;
  /* a low pathogenic site: one basic residue before HA2 */
  cleavage_counts_clear(&c);
  add_n(&c, LPAI_NT, 50);
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_str(r.aa, "NVPARSSR");
  test_eq_size(r.reading, CLEAVAGE_MONOBASIC);
  test_eq_size(r.verdict, CLEAVAGE_FIRM);
  /* a motif with a stop codon cannot be read */
  cleavage_counts_clear(&c);
  add_n(&c, "CAAAGAGAGTAGAGAAGAAAAAGA", 50);
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_size(r.reading, CLEAVAGE_UNREADABLE);
  /* and one with an ambiguous base */
  cleavage_counts_clear(&c);
  add_n(&c, "CAAAGAGAGAGAAGAAGAAAAANA", 50);
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_size(r.reading, CLEAVAGE_UNREADABLE);
}

static void test_overflow(void) {
  /* more distinct motifs than the table holds: the reads are counted,
     the support is a lower bound, and the frequent motif is kept */
  cleavage_counts_t c;
  cleavage_counts_clear(&c);
  add_n(&c, HPAI_NT, 500);
  /* the other motifs are built on a different template, so that none of
     them is the frequent one */
  char nt[CLEAVAGE_SITE_NT + 1];
  memcpy(nt, LPAI_NT, sizeof(nt));
  for (size_t i = 0; i < CLEAVAGE_MAX_MOTIFS + 100; i++) {
    for (size_t j = 0; j < 6; j++)
      nt[j] = "ACGT"[(i >> (2 * j)) & 3];
    add_n(&c, nt, 1);
  }
  test_eq_size(c.n_motifs, CLEAVAGE_MAX_MOTIFS);
  test_assert(c.overflow > 0);
  cleavage_result_t r;
  cleavage_result_of(&c, CLEAVAGE_CHEMISTRY_R10, &r);
  test_eq_str(r.nt, HPAI_NT);
  test_eq_size(r.motif_reads, 500);
}

static void test_save_load(void) {
  cleavage_counts_t c;
  cleavage_counts_clear(&c);
  add_n(&c, HPAI_NT, 30);
  add_n(&c, LPAI_NT, 7);
  c.overflow = 3;
  char *buf = NULL;
  size_t size = 0;
  FILE *f = open_memstream(&buf, &size);
  test_assert(f != NULL);
  cleavage_counts_save(&c, f);
  fclose(f);
  cleavage_counts_t loaded;
  cleavage_counts_clear(&loaded);
  for (char *line = strtok(buf, "\n"); line != NULL;
       line = strtok(NULL, "\n"))
    test_assert(cleavage_counts_load_line(&loaded, line));
  test_eq_size(loaded.searched, 37);
  test_eq_size(loaded.overflow, 3);
  test_eq_size(cleavage_counts_reads(&loaded), 40);
  test_eq_size(loaded.n_motifs, 2);
  free(buf);
  /* a line of another kind is not one of these */
  char other[] = "reads 1 2 3";
  test_assert(!cleavage_counts_load_line(&loaded, other));
}

static void test_chemistry(void) {
  test_eq_size(cleavage_chemistry(true, NULL), CLEAVAGE_CHEMISTRY_ILLUMINA);
  test_eq_size(cleavage_chemistry(true, "dna_r9.4.1_e8_hac@v3.3"),
               CLEAVAGE_CHEMISTRY_ILLUMINA);
  test_eq_size(cleavage_chemistry(false, NULL), CLEAVAGE_CHEMISTRY_UNKNOWN);
  test_eq_size(cleavage_chemistry(false, ""), CLEAVAGE_CHEMISTRY_UNKNOWN);
  test_eq_size(cleavage_chemistry(false, "dna_r10.4.1_e8.2_400bps_hac@v4.2.0"),
               CLEAVAGE_CHEMISTRY_R10);
  test_eq_size(cleavage_chemistry(false, "dna_r9.4.1_e8_hac@v3.3"),
               CLEAVAGE_CHEMISTRY_R9);
  test_eq_size(cleavage_chemistry(false, "something_else"),
               CLEAVAGE_CHEMISTRY_UNKNOWN);
  test_eq_str(cleavage_chemistry_name(CLEAVAGE_CHEMISTRY_R10), "R10");
  test_eq_str(cleavage_reading_name(CLEAVAGE_MULTIBASIC), "multibasic");
  test_eq_str(cleavage_verdict_name(CLEAVAGE_CONFIRM), "confirm");
}

void cleavage_site_test(void) {
  test_counts();
  test_supported_readings();
  test_control();
  test_fusion();
  test_readings();
  test_overflow();
  test_save_load();
  test_chemistry();
}
