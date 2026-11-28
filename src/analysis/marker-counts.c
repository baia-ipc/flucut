/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "marker-counts.h"

/* --- the tally of one window --- */

/* FNV-1a over a value */
static size_t hash(const char *value) {
  size_t h = 2166136261u;
  for (const char *p = value; *p != '\0'; p++) {
    h ^= (unsigned char)*p;
    h *= 16777619u;
  }
  return h;
}

/* the slot of <value>: the one holding it, or the first empty one after
   its place */
static size_t slot_of(const marker_tally_t *t, const char *value) {
  size_t s = hash(value) & (MARKER_SLOTS - 1);
  while (t->slots[s] != 0 &&
         strcmp(t->entries[t->slots[s] - 1].value, value) != 0)
    s = (s + 1) & (MARKER_SLOTS - 1);
  return s;
}

static void tally_add(marker_tally_t *t, const char *value, size_t reads) {
  if (value[0] == '\0' || strlen(value) > MARKER_VALUE_MAX)
    return;
  t->reads += reads;
  size_t s = slot_of(t, value);
  if (t->slots[s] != 0) {
    t->entries[t->slots[s] - 1].reads += reads;
    return;
  }
  if (t->n == MARKER_VALUES) {
    t->overflow += reads;
    return;
  }
  snprintf(t->entries[t->n].value, sizeof(t->entries[t->n].value), "%s",
           value);
  t->entries[t->n].reads = reads;
  t->n++;
  t->slots[s] = t->n;
}

static void tally_merge(marker_tally_t *t, const marker_tally_t *other) {
  t->overflow += other->overflow;
  t->reads += other->overflow;
  for (size_t i = 0; i < other->n; i++)
    tally_add(t, other->entries[i].value, other->entries[i].reads);
}

/* --- the counts --- */

void marker_counts_clear(marker_counts_t *self) {
  assert(self != NULL);
  memset(self, 0, sizeof(marker_counts_t));
}

void marker_counts_set_reference(marker_counts_t *self, size_t marker,
    const char *value, const char *anchor_aa) {
  assert(self != NULL);
  assert(marker < MARKER_N_MARKERS);
  assert(value != NULL);
  assert(anchor_aa != NULL);
  marker_one_counts_t *m = self->markers + marker;
  if (m->reference_value[0] != '\0')
    return;
  snprintf(m->reference_value, sizeof(m->reference_value), "%s", value);
  snprintf(m->reference_anchor, sizeof(m->reference_anchor), "%s", anchor_aa);
}

void marker_counts_searched(marker_counts_t *self, size_t marker) {
  assert(self != NULL);
  assert(marker < MARKER_N_MARKERS);
  self->markers[marker].searched++;
}

void marker_counts_add(marker_counts_t *self, size_t marker,
    const char *value) {
  assert(self != NULL);
  assert(marker < MARKER_N_MARKERS);
  assert(value != NULL);
  tally_add(&self->markers[marker].value, value, 1);
}

void marker_counts_add_anchor(marker_counts_t *self, size_t marker,
    const char *aa) {
  assert(self != NULL);
  assert(marker < MARKER_N_MARKERS);
  assert(aa != NULL);
  tally_add(&self->markers[marker].anchor, aa, 1);
}

size_t marker_counts_reads(const marker_counts_t *self, size_t marker) {
  assert(self != NULL);
  assert(marker < MARKER_N_MARKERS);
  return self->markers[marker].value.reads;
}

bool marker_counts_any(const marker_counts_t *self) {
  assert(self != NULL);
  for (size_t m = 0; m < MARKER_N_MARKERS; m++)
    if (self->markers[m].searched > 0)
      return true;
  return false;
}

void marker_counts_merge(marker_counts_t *self, const marker_counts_t *other) {
  assert(self != NULL);
  assert(other != NULL);
  for (size_t m = 0; m < MARKER_N_MARKERS; m++) {
    const marker_one_counts_t *o = other->markers + m;
    self->markers[m].searched += o->searched;
    tally_merge(&self->markers[m].value, &o->value);
    tally_merge(&self->markers[m].anchor, &o->anchor);
    if (o->reference_value[0] != '\0')
      marker_counts_set_reference(self, m, o->reference_value,
                                  o->reference_anchor);
  }
}

/* --- saved and loaded --- */

void marker_counts_save(const marker_counts_t *self, FILE *file) {
  assert(self != NULL);
  assert(file != NULL);
  for (size_t m = 0; m < MARKER_N_MARKERS; m++) {
    const marker_one_counts_t *c = self->markers + m;
    if (c->searched == 0)
      continue;
    fprintf(file, "marker-reads %s %zu %zu %zu\n", marker_def(m)->name,
            c->searched, c->value.overflow, c->anchor.overflow);
    if (c->reference_value[0] != '\0')
      fprintf(file, "marker-reference %s %s %s\n", marker_def(m)->name,
              c->reference_value, c->reference_anchor);
    for (size_t i = 0; i < c->value.n; i++)
      fprintf(file, "marker-value %s %s %zu\n", marker_def(m)->name,
              c->value.entries[i].value, c->value.entries[i].reads);
    for (size_t i = 0; i < c->anchor.n; i++)
      fprintf(file, "marker-anchor %s %s %zu\n", marker_def(m)->name,
              c->anchor.entries[i].value, c->anchor.entries[i].reads);
  }
}

/* the marker named at the start of <rest>, and <rest> moved past it;
   MARKER_N_MARKERS if it names none */
static size_t marker_named(char **rest) {
  char name[32];
  int used = 0;
  if (sscanf(*rest, "%31s%n", name, &used) != 1)
    return MARKER_N_MARKERS;
  for (size_t m = 0; m < MARKER_N_MARKERS; m++)
    if (strcmp(marker_def(m)->name, name) == 0) {
      *rest += used;
      return m;
    }
  return MARKER_N_MARKERS;
}

bool marker_counts_load_line(marker_counts_t *self, char *line) {
  assert(self != NULL);
  assert(line != NULL);
  char *rest;
  if (strncmp(line, "marker-reads ", 13) == 0) {
    rest = line + 13;
    size_t m = marker_named(&rest);
    if (m == MARKER_N_MARKERS)
      return false;
    char *next;
    self->markers[m].searched += strtoull(rest, &next, 10);
    size_t over = strtoull(next, &next, 10);
    self->markers[m].value.overflow += over;
    self->markers[m].value.reads += over;
    over = strtoull(next, NULL, 10);
    self->markers[m].anchor.overflow += over;
    self->markers[m].anchor.reads += over;
    return true;
  }
  if (strncmp(line, "marker-reference ", 17) == 0) {
    rest = line + 17;
    size_t m = marker_named(&rest);
    char value[MARKER_VALUE_MAX + 2], anchor[MARKER_ANCHOR_AA + 2] = "";
    if (m == MARKER_N_MARKERS || sscanf(rest, "%21s %9s", value, anchor) < 1 ||
        strlen(value) > MARKER_VALUE_MAX)
      return false;
    marker_counts_set_reference(self, m, value, anchor);
    return true;
  }
  bool is_value = strncmp(line, "marker-value ", 13) == 0;
  if (!is_value && strncmp(line, "marker-anchor ", 14) != 0)
    return false;
  rest = line + (is_value ? 13 : 14);
  size_t m = marker_named(&rest);
  char value[MARKER_VALUE_MAX + 2];
  int used = 0;
  if (m == MARKER_N_MARKERS || sscanf(rest, "%21s%n", value, &used) != 1 ||
      strlen(value) > MARKER_VALUE_MAX)
    return false;
  size_t reads = strtoull(rest + used, NULL, 10);
  tally_add(is_value ? &self->markers[m].value : &self->markers[m].anchor,
            value, reads);
  return true;
}

/* --- what the counts say --- */

const char* marker_verdict_name(size_t verdict) {
  switch (verdict) {
    case MARKER_FIRM: return "firm";
    case MARKER_CONFIRM: return "confirm";
    case MARKER_TOO_FEW_READS: return "too few reads";
    default: return "unknown";
  }
}

/* the entry of <t> with the most reads, or -1 if it has none */
static long most_frequent(const marker_tally_t *t) {
  long best = -1;
  for (size_t i = 0; i < t->n; i++)
    if (best < 0 || t->entries[i].reads > t->entries[best].reads)
      best = (long)i;
  return best;
}

/* the share of the reads of <t> agreeing on its most frequent value (0
   where it has no reads) */
static double agreement(const marker_tally_t *t) {
  long best = most_frequent(t);
  if (best < 0 || t->reads == 0)
    return 0;
  return (double)t->entries[best].reads / (double)t->reads;
}

void marker_result_of(const marker_counts_t *counts, size_t marker,
    marker_result_t *out) {
  assert(counts != NULL);
  assert(out != NULL);
  assert(marker < MARKER_N_MARKERS);
  const marker_def_t *def = marker_def(marker);
  const marker_one_counts_t *c = counts->markers + marker;
  memset(out, 0, sizeof(marker_result_t));
  out->searched = c->searched;
  out->reads = c->value.reads;
  out->distinct = c->value.n;
  out->anchor_reads = c->anchor.reads;
  out->anchor_support = agreement(&c->anchor);
  out->anchor_weak = out->anchor_reads > 0 &&
                     out->anchor_support < MARKER_MIN_ANCHOR_SUPPORT;
  snprintf(out->reference_value, sizeof(out->reference_value), "%s",
           c->reference_value);
  long best = most_frequent(&c->value);
  if (best < 0 || out->reads == 0) {
    out->verdict = MARKER_TOO_FEW_READS;
    return;
  }
  out->found = true;
  snprintf(out->value, sizeof(out->value), "%s", c->value.entries[best].value);
  out->value_reads = c->value.entries[best].reads;
  out->support = (double)out->value_reads / (double)out->reads;
  out->differs_from_reference = out->reference_value[0] != '\0' &&
      strcmp(out->value, out->reference_value) != 0;
  /* the values the sample holds, and those of them worth looking at */
  for (size_t i = 0; i < c->value.n; i++) {
    size_t reads = c->value.entries[i].reads;
    double share = (double)reads / (double)out->reads;
    if (share < MARKER_MIN_SHARE || reads < MARKER_MIN_VALUE_READS)
      continue;
    out->supported++;
    if (marker_value_attention(def, c->value.entries[i].value)) {
      out->attention = true;
      out->attention_share += share;
    }
  }
  /* the most frequent values, most frequent first */
  bool taken[MARKER_VALUES] = {false};
  for (size_t k = 0; k < MARKER_TOP_VALUES; k++) {
    long at = -1;
    for (size_t i = 0; i < c->value.n; i++)
      if (!taken[i] && (at < 0 ||
                        c->value.entries[i].reads > c->value.entries[at].reads))
        at = (long)i;
    if (at < 0 || c->value.entries[at].reads == 0)
      break;
    taken[at] = true;
    marker_value_t *v = out->top + out->n_top++;
    snprintf(v->value, sizeof(v->value), "%s", c->value.entries[at].value);
    v->reads = c->value.entries[at].reads;
    v->share = (double)v->reads / (double)out->reads;
    v->attention = marker_value_attention(def, v->value);
  }
  out->verdict = out->reads < MARKER_MIN_READS ? MARKER_TOO_FEW_READS :
      (out->support >= MARKER_MIN_SUPPORT && !out->anchor_weak ?
       MARKER_FIRM : MARKER_CONFIRM);
}
