/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "counter.h"
#include "logger.h"
#include "test.h"

struct counter_t {
  counter_reads_t reads;
  /* reads by primer found, for each end */
  size_t n_variants[COUNTER_N_ENDS];
  size_t *variants[COUNTER_N_ENDS];
  /* internal ends: for each reference, NULL or the counts by end and
     position (end * length + position) */
  const fasta_t *refs;
  size_t **sites;
  lengths_t *lengths;
};

size_t counter_n_length_groups(const fasta_t *references) {
  return references == NULL ? 1 : references->n_records + 1;
}

counter_t* counter_new(const scheme_t *scheme, const fasta_t *references) {
  assert(scheme != NULL);
  counter_t *self;
  XCALLOC(self, 1, sizeof(counter_t));
  self->n_variants[COUNTER_LEFT] = scheme->n_left;
  self->n_variants[COUNTER_RIGHT] = scheme->n_right;
  for (size_t e = 0; e < COUNTER_N_ENDS; e++)
    XCALLOC(self->variants[e], self->n_variants[e] + 1, sizeof(size_t));
  self->refs = references;
  if (references != NULL)
    XCALLOC(self->sites, references->n_records + 1, sizeof(size_t *));
  self->lengths = lengths_new(counter_n_length_groups(references));
  return self;
}

void counter_delete(counter_t *self) {
  if (self == NULL)
    return;
  for (size_t e = 0; e < COUNTER_N_ENDS; e++)
    free(self->variants[e]);
  if (self->refs != NULL)
    for (size_t r = 0; r < self->refs->n_records; r++)
      free(self->sites[r]);
  free(self->sites);
  lengths_delete(self->lengths);
  free(self);
}

void counter_clear(counter_t *self) {
  assert(self != NULL);
  memset(&self->reads, 0, sizeof(counter_reads_t));
  for (size_t e = 0; e < COUNTER_N_ENDS; e++)
    memset(self->variants[e], 0, sizeof(size_t) * self->n_variants[e]);
  if (self->refs != NULL)
    for (size_t r = 0; r < self->refs->n_records; r++)
      if (self->sites[r] != NULL)
        memset(self->sites[r], 0, sizeof(size_t) * COUNTER_N_ENDS *
               self->refs->records[r].seq_len);
  lengths_clear(self->lengths);
}

/* the counts of the sites of the reference <r> (allocated at its first
   site) */
static size_t* ref_sites(counter_t *self, size_t r) {
  if (self->sites[r] == NULL)
    XCALLOC(self->sites[r], COUNTER_N_ENDS * self->refs->records[r].seq_len,
            sizeof(size_t));
  return self->sites[r];
}

void counter_merge(counter_t *self, const counter_t *other) {
  assert(self != NULL);
  assert(other != NULL);
  assert(self->refs == other->refs);
  counter_reads_t *s = &self->reads;
  const counter_reads_t *o = &other->reads;
  s->reads_in += o->reads_in;
  s->bases_in += o->bases_in;
  s->reads_out += o->reads_out;
  s->bases_out += o->bases_out;
  s->too_short += o->too_short;
  s->dimers += o->dimers;
  for (size_t i = 0; i < 3; i++)
    s->n_ends[i] += o->n_ends[i];
  s->minus += o->minus;
  s->with_internal += o->with_internal;
  s->adapter_start += o->adapter_start;
  s->adapter_end += o->adapter_end;
  for (size_t e = 0; e < COUNTER_N_ENDS; e++) {
    for (size_t i = 0; i < 3; i++)
      s->origin[e][i] += o->origin[e][i];
    assert(self->n_variants[e] == other->n_variants[e]);
    for (size_t v = 0; v < self->n_variants[e]; v++)
      self->variants[e][v] += other->variants[e][v];
  }
  if (self->refs != NULL)
    for (size_t r = 0; r < self->refs->n_records; r++) {
      if (other->sites[r] == NULL)
        continue;
      size_t *dest = ref_sites(self, r);
      size_t n = COUNTER_N_ENDS * self->refs->records[r].seq_len;
      for (size_t i = 0; i < n; i++)
        dest[i] += other->sites[r][i];
    }
  lengths_add_all(self->lengths, other->lengths);
}

/* count an end found (<e>: COUNTER_LEFT or COUNTER_RIGHT) */
static void count_end(counter_t *self, const amplicon_end_t *end,
    size_t e) {
  if (!end->found)
    return;
  self->reads.origin[e][end->origin]++;
  assert(end->variant < self->n_variants[e]);
  self->variants[e][end->variant]++;
  if (end->origin != ORIGIN_INTERNAL || !end->ref.mapped ||
      self->refs == NULL)
    return;
  size_t r = end->ref.ref_index;
  assert(r < self->refs->n_records);
  size_t len = self->refs->records[r].seq_len;
  long pos = end->ref.ref_start;
  if (pos < 0 || (size_t)pos >= len)
    return;
  ref_sites(self, r)[e * len + (size_t)pos]++;
}

/* the counts of counter_reads_t, in the order of the files of counts:
   10 single counts, n_ends, origin */
#define N_READ_COUNTS (10 + 3 + COUNTER_N_ENDS * 3)

static void read_counts(counter_reads_t *r, size_t *out[N_READ_COUNTS]) {
  size_t *single[10] = {&r->reads_in, &r->bases_in, &r->reads_out,
      &r->bases_out, &r->too_short, &r->dimers, &r->minus,
      &r->with_internal, &r->adapter_start, &r->adapter_end};
  size_t n = 0;
  for (size_t i = 0; i < 10; i++)
    out[n++] = single[i];
  for (size_t i = 0; i < 3; i++)
    out[n++] = r->n_ends + i;
  for (size_t e = 0; e < COUNTER_N_ENDS; e++)
    for (size_t i = 0; i < 3; i++)
      out[n++] = r->origin[e] + i;
  assert(n == N_READ_COUNTS);
}

#define COUNTS_HEADER "flucut-counts 1"

void counter_save(const counter_t *self, FILE *file) {
  assert(self != NULL);
  assert(file != NULL);
  fputs(COUNTS_HEADER "\nreads", file);
  counter_reads_t copy = self->reads;
  size_t *p[N_READ_COUNTS];
  read_counts(&copy, p);
  for (size_t i = 0; i < N_READ_COUNTS; i++)
    fprintf(file, " %zu", *p[i]);
  fputc('\n', file);
  for (size_t e = 0; e < COUNTER_N_ENDS; e++) {
    fprintf(file, "variants %zu %zu", e, self->n_variants[e]);
    for (size_t v = 0; v < self->n_variants[e]; v++)
      fprintf(file, " %zu", self->variants[e][v]);
    fputc('\n', file);
  }
  if (self->refs != NULL)
    for (size_t r = 0; r < self->refs->n_records; r++) {
      if (self->sites[r] == NULL)
        continue;
      size_t len = self->refs->records[r].seq_len;
      for (size_t i = 0; i < COUNTER_N_ENDS * len; i++)
        if (self->sites[r][i] > 0)
          fprintf(file, "site %zu %zu %zu %zu\n", r, i / len, i % len,
                  self->sites[r][i]);
    }
  for (size_t g = 0; g < lengths_n_groups(self->lengths); g++)
    for (size_t b = 0; b < LENGTHS_N_BINS; b++) {
      size_t n = lengths_count(self->lengths, g, b);
      if (n > 0)
        fprintf(file, "length %zu %zu %zu\n", g, b, n);
    }
  fputs("end\n", file);
}

/* load a line of a file of counts (see counter_load()) */
static bool load_line(counter_t *self, char *line, bool *ended) {
  char *rest;
  if (strncmp(line, "reads", 5) == 0) {
    size_t *p[N_READ_COUNTS];
    read_counts(&self->reads, p);
    rest = line + 5;
    for (size_t i = 0; i < N_READ_COUNTS; i++) {
      char *end;
      unsigned long long v = strtoull(rest, &end, 10);
      if (end == rest)
        return false;
      *p[i] += (size_t)v;
      rest = end;
    }
    return true;
  }
  size_t a, b, c, d;
  int n;
  if (sscanf(line, "variants %zu %zu%n", &a, &b, &n) == 2) {
    if (a >= COUNTER_N_ENDS || b != self->n_variants[a])
      return false;
    rest = line + n;
    for (size_t v = 0; v < b; v++) {
      char *end;
      unsigned long long x = strtoull(rest, &end, 10);
      if (end == rest)
        return false;
      self->variants[a][v] += (size_t)x;
      rest = end;
    }
    return true;
  }
  if (sscanf(line, "site %zu %zu %zu %zu", &a, &b, &c, &d) == 4) {
    if (self->refs == NULL || a >= self->refs->n_records ||
        b >= COUNTER_N_ENDS || c >= self->refs->records[a].seq_len)
      return false;
    ref_sites(self, a)[b * self->refs->records[a].seq_len + c] += d;
    return true;
  }
  if (sscanf(line, "length %zu %zu %zu", &a, &b, &c) == 3) {
    if (a >= lengths_n_groups(self->lengths) || b >= LENGTHS_N_BINS)
      return false;
    lengths_add_to_bin(self->lengths, a, b, c);
    return true;
  }
  if (strcmp(line, "end\n") == 0) {
    *ended = true;
    return true;
  }
  return false;
}

bool counter_load(counter_t *self, FILE *file) {
  assert(self != NULL);
  assert(file != NULL);
  char *line = NULL;
  size_t cap = 0;
  bool ok = getline(&line, &cap, file) > 0 &&
            strcmp(line, COUNTS_HEADER "\n") == 0, ended = false;
  while (ok && !ended && getline(&line, &cap, file) > 0)
    ok = load_line(self, line, &ended);
  free(line);
  return ok && ended;
}

/* the length of the template of a complete amplicon, in its group (see
   counter_lengths()) */
static void count_length(counter_t *self, const amplicon_t *a) {
  if (!a->left.found || !a->right.found)
    return;
  long len = a->right.hit.start - a->left.hit.end - 1;
  if (len <= 0)
    return;
  size_t group = 0;
  if (self->refs != NULL) {
    group = self->refs->n_records;
    if (a->left.ref.mapped)
      group = a->left.ref.ref_index;
    else if (a->right.ref.mapped)
      group = a->right.ref.ref_index;
  }
  lengths_add(self->lengths, group, (size_t)len);
}

void counter_add_read(counter_t *self, const amplicon_t *a, size_t len) {
  assert(self != NULL);
  assert(a != NULL);
  counter_reads_t *s = &self->reads;
  s->reads_in++;
  s->bases_in += len;
  size_t n_ends = amplicon_n_ends(a);
  s->n_ends[n_ends]++;
  if (a->strand == STRAND_MINUS)
    s->minus++;
  if (n_ends == 2 && amplicon_kept_length(a) == 0)
    s->dimers++;
  if (a->n_internal > 0)
    s->with_internal++;
  if (a->adapter_start.found)
    s->adapter_start++;
  if (a->adapter_end.found)
    s->adapter_end++;
  count_end(self, &a->left, COUNTER_LEFT);
  count_end(self, &a->right, COUNTER_RIGHT);
  count_length(self, a);
}

void counter_add_written(counter_t *self, size_t len) {
  self->reads.reads_out++;
  self->reads.bases_out += len;
}

void counter_add_too_short(counter_t *self) {
  self->reads.too_short++;
}

const counter_reads_t* counter_reads(const counter_t *self) {
  assert(self != NULL);
  return &self->reads;
}

size_t counter_n_variants(const counter_t *self, size_t end) {
  assert(end < COUNTER_N_ENDS);
  return self->n_variants[end];
}

size_t counter_primer(const counter_t *self, size_t end, size_t variant) {
  assert(end < COUNTER_N_ENDS);
  assert(variant < self->n_variants[end]);
  return self->variants[end][variant];
}

const size_t* counter_sites(const counter_t *self, size_t ref, size_t end) {
  assert(end < COUNTER_N_ENDS);
  if (self->refs == NULL)
    return NULL;
  assert(ref < self->refs->n_records);
  if (self->sites[ref] == NULL)
    return NULL;
  return self->sites[ref] + end * self->refs->records[ref].seq_len;
}

const lengths_t* counter_lengths(const counter_t *self) {
  return self->lengths;
}

/* --- unit tests --- */

/* an amplicon with both ends found, from <from> to <to> (the template is
   between them), the left end internal at <site> of reference 1 */
static amplicon_t test_amplicon(long from, long to, long site) {
  amplicon_t a;
  memset(&a, 0, sizeof(a));
  a.strand = STRAND_PLUS;
  a.left.found = a.right.found = true;
  a.left.variant = 1;
  a.left.hit.end = from - 1;
  a.right.hit.start = to + 1;
  a.left.origin = ORIGIN_INTERNAL;
  a.left.ref.mapped = true;
  a.left.ref.ref_index = 1;
  a.left.ref.ref_start = site;
  a.right.origin = ORIGIN_TERMINAL;
  a.keep_from = from;
  a.keep_to = to;
  return a;
}

void counter_test(void) {
  Logger *logger = new_logger(stderr, LOG_ERROR);
  scheme_t *sc = scheme_new_builtin("mbt", logger);
  fasta_record_t recs[2] = {{"r0", "", "ACGTACGTAC", 10},
                            {"r1", "", "ACGTACGTACGTACGTACGT", 20}};
  fasta_t refs = {recs, 2};
  counter_t *c = counter_new(sc, &refs), *sum = counter_new(sc, &refs);
  test_eq_size(counter_n_variants(c, COUNTER_LEFT), sc->n_left);
  test_eq_size(counter_n_variants(c, COUNTER_RIGHT), sc->n_right);
  test_eq_size(counter_n_length_groups(&refs), 3);
  test_eq_size(counter_n_length_groups(NULL), 1);
  test_assert(counter_sites(c, 1, COUNTER_LEFT) == NULL);
  amplicon_t a = test_amplicon(30, 129, 7);
  counter_add_read(c, &a, 200);
  counter_add_written(c, 100);
  a = test_amplicon(30, 29, 7);
  counter_add_read(c, &a, 60);
  counter_add_too_short(c);
  const counter_reads_t *r = counter_reads(c);
  test_eq_size(r->reads_in, 2);
  test_eq_size(r->bases_in, 260);
  test_eq_size(r->reads_out, 1);
  test_eq_size(r->n_ends[2], 2);
  test_eq_size(r->dimers, 1);
  test_eq_size(r->too_short, 1);
  test_eq_size(r->origin[COUNTER_LEFT][ORIGIN_INTERNAL], 2);
  test_eq_size(r->origin[COUNTER_RIGHT][ORIGIN_TERMINAL], 2);
  test_eq_size(counter_primer(c, COUNTER_LEFT, 1), 2);
  test_eq_size(counter_primer(c, COUNTER_LEFT, 0), 0);
  test_eq_size(counter_primer(c, COUNTER_RIGHT, 0), 2);
  const size_t *s = counter_sites(c, 1, COUNTER_LEFT);
  test_assert(s != NULL);
  test_eq_size(s[7], 2);
  test_eq_size(s[6], 0);
  test_eq_size(counter_sites(c, 1, COUNTER_RIGHT)[7], 0);
  test_assert(counter_sites(c, 0, COUNTER_LEFT) == NULL);
  /* the length of the template (100), in the group of reference 1; the
     dimer has none */
  test_eq_size(lengths_count(counter_lengths(c), 1, 10), 1);
  /* the sum of two counters */
  counter_merge(sum, c);
  counter_merge(sum, c);
  test_eq_size(counter_reads(sum)->reads_in, 4);
  test_eq_size(counter_primer(sum, COUNTER_LEFT, 1), 4);
  test_eq_size(counter_sites(sum, 1, COUNTER_LEFT)[7], 4);
  test_eq_size(lengths_count(counter_lengths(sum), 1, 10), 2);
  /* saved and loaded (added): the same counts, twice */
  FILE *f = tmpfile();
  test_assert(f != NULL);
  counter_save(sum, f);
  rewind(f);
  counter_t *loaded = counter_new(sc, &refs);
  test_assert(counter_load(loaded, f));
  rewind(f);
  test_assert(counter_load(loaded, f));
  test_eq_size(counter_reads(loaded)->reads_in, 8);
  test_eq_size(counter_reads(loaded)->origin[COUNTER_LEFT][ORIGIN_INTERNAL],
               8);
  test_eq_size(counter_primer(loaded, COUNTER_LEFT, 1), 8);
  test_eq_size(counter_sites(loaded, 1, COUNTER_LEFT)[7], 8);
  test_eq_size(lengths_count(counter_lengths(loaded), 1, 10), 4);
  /* not a file of counts, or incomplete */
  rewind(f);
  fputs("something else\n", f);
  rewind(f);
  test_assert(!counter_load(loaded, f));
  fclose(f);
  f = tmpfile();
  fputs("flucut-counts 1\nreads 1 2\n", f);
  rewind(f);
  test_assert(!counter_load(loaded, f));
  fclose(f);
  counter_delete(loaded);
  counter_clear(sum);
  test_eq_size(counter_reads(sum)->reads_in, 0);
  test_eq_size(counter_sites(sum, 1, COUNTER_LEFT)[7], 0);
  test_eq_size(counter_primer(sum, COUNTER_LEFT, 1), 0);
  counter_delete(c);
  counter_delete(sum);
  scheme_delete(sc);
  delete_logger(logger);
}
