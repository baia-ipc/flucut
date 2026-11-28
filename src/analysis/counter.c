/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "counter.h"
#include "marker-read.h"
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
  /* reconstructed coverage: for each reference, NULL or an array of
     2 * COVER_N_CATEGORIES * (seq_len + 1) counts; the start counts
     (side 0) and the end counts (side 1, the position after the span) of
     each category, prefix-summed by counter_coverage() */
  size_t **cover;
  lengths_t *lengths;
  read_stats_t *stats;
  read_metadata_t metadata;
  /* the cleavage site of the haemagglutinin read from the reads */
  cleavage_counts_t cleavage;
  /* the markers read from the reads (markers.h) */
  marker_counts_t markers;
  /* the reads of each HA and NA subtype among the first reads */
  subtype_reads_t subtype_reads;
};

/* the index in cover[r] of the start (side 0) or end (side 1) count of
   category <cat> at position <pos>, for a reference of length <len> */
#define COVER_SLOTS(len) (COVER_N_CATEGORIES * ((len) + 1))
#define COVER_IDX(side, cat, pos, len) \
  ((side) * COVER_SLOTS(len) + (cat) * ((len) + 1) + (pos))

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
  if (references != NULL) {
    XCALLOC(self->sites, references->n_records + 1, sizeof(size_t *));
    XCALLOC(self->cover, references->n_records + 1, sizeof(size_t *));
  }
  self->lengths = lengths_new(counter_n_length_groups(references));
  XCALLOC(self->stats, 1, sizeof(read_stats_t));
  static const char *const stats_names[COUNTER_STATS_N_SETS] = {
    "input", "written"};
  read_stats_init(self->stats, COUNTER_STATS_N_SETS, stats_names);
  return self;
}

void counter_delete(counter_t *self) {
  if (self == NULL)
    return;
  for (size_t e = 0; e < COUNTER_N_ENDS; e++)
    free(self->variants[e]);
  if (self->refs != NULL)
    for (size_t r = 0; r < self->refs->n_records; r++) {
      free(self->sites[r]);
      free(self->cover[r]);
    }
  free(self->sites);
  free(self->cover);
  lengths_delete(self->lengths);
  read_stats_free(self->stats);
  free(self->stats);
  free(self);
}

void counter_clear(counter_t *self) {
  assert(self != NULL);
  memset(&self->reads, 0, sizeof(counter_reads_t));
  for (size_t e = 0; e < COUNTER_N_ENDS; e++)
    memset(self->variants[e], 0, sizeof(size_t) * self->n_variants[e]);
  if (self->refs != NULL)
    for (size_t r = 0; r < self->refs->n_records; r++) {
      if (self->sites[r] != NULL)
        memset(self->sites[r], 0, sizeof(size_t) * COUNTER_N_ENDS *
               self->refs->records[r].seq_len);
      if (self->cover[r] != NULL)
        memset(self->cover[r], 0, sizeof(size_t) * 2 *
               COVER_SLOTS(self->refs->records[r].seq_len));
    }
  lengths_clear(self->lengths);
  read_stats_clear(self->stats);
  read_metadata_clear(&self->metadata);
  cleavage_counts_clear(&self->cleavage);
  marker_counts_clear(&self->markers);
  subtype_reads_clear(&self->subtype_reads);
}

/* the counts of the sites of the reference <r> (allocated at its first
   site) */
static size_t* ref_sites(counter_t *self, size_t r) {
  if (self->sites[r] == NULL)
    XCALLOC(self->sites[r], COUNTER_N_ENDS * self->refs->records[r].seq_len,
            sizeof(size_t));
  return self->sites[r];
}

/* the coverage counts of the reference <r> (allocated at its first read) */
static size_t* ref_cover(counter_t *self, size_t r) {
  if (self->cover[r] == NULL)
    XCALLOC(self->cover[r], 2 * COVER_SLOTS(self->refs->records[r].seq_len),
            sizeof(size_t));
  return self->cover[r];
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
  for (size_t k = 0; k < READ_N_CATEGORIES; k++) {
    s->written[k] += o->written[k];
    s->filtered[k] += o->filtered[k];
  }
  s->too_short += o->too_short;
  s->dimers += o->dimers;
  for (size_t i = 0; i < 3; i++)
    s->n_ends[i] += o->n_ends[i];
  s->minus += o->minus;
  s->with_internal += o->with_internal;
  s->split += o->split;
  s->cut += o->cut;
  s->cut_again += o->cut_again;
  s->adapter_start += o->adapter_start;
  s->adapter_end += o->adapter_end;
  s->placed += o->placed;
  for (size_t e = 0; e < COUNTER_N_ENDS; e++) {
    for (size_t i = 0; i < 3; i++)
      s->origin[e][i] += o->origin[e][i];
    assert(self->n_variants[e] == other->n_variants[e]);
    for (size_t v = 0; v < self->n_variants[e]; v++)
      self->variants[e][v] += other->variants[e][v];
  }
  cleavage_counts_merge(&self->cleavage, &other->cleavage);
  marker_counts_merge(&self->markers, &other->markers);
  subtype_reads_merge(&self->subtype_reads, &other->subtype_reads);
  if (self->refs != NULL)
    for (size_t r = 0; r < self->refs->n_records; r++) {
      if (other->sites[r] == NULL)
        continue;
      size_t *dest = ref_sites(self, r);
      size_t n = COUNTER_N_ENDS * self->refs->records[r].seq_len;
      for (size_t i = 0; i < n; i++)
        dest[i] += other->sites[r][i];
    }
  if (self->refs != NULL)
    for (size_t r = 0; r < self->refs->n_records; r++) {
      if (other->cover[r] == NULL)
        continue;
      size_t *dest = ref_cover(self, r);
      size_t n = 2 * COVER_SLOTS(self->refs->records[r].seq_len);
      for (size_t i = 0; i < n; i++)
        dest[i] += other->cover[r][i];
    }
  lengths_add_all(self->lengths, other->lengths);
  read_stats_merge(self->stats, other->stats);
  read_metadata_merge(&self->metadata, &other->metadata);
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
   14 single counts, n_ends, origin, written, filtered */
#define N_SINGLE_COUNTS 14
#define N_READ_COUNTS (N_SINGLE_COUNTS + 3 + COUNTER_N_ENDS * 3 + \
                       2 * READ_N_CATEGORIES)

static void read_counts(counter_reads_t *r, size_t *out[N_READ_COUNTS]) {
  size_t *single[N_SINGLE_COUNTS] = {&r->reads_in, &r->bases_in,
      &r->reads_out, &r->bases_out, &r->too_short, &r->dimers, &r->minus,
      &r->with_internal, &r->split, &r->adapter_start, &r->adapter_end,
      &r->cut, &r->placed, &r->cut_again};
  size_t n = 0;
  for (size_t i = 0; i < N_SINGLE_COUNTS; i++)
    out[n++] = single[i];
  for (size_t i = 0; i < 3; i++)
    out[n++] = r->n_ends + i;
  for (size_t e = 0; e < COUNTER_N_ENDS; e++)
    for (size_t i = 0; i < 3; i++)
      out[n++] = r->origin[e] + i;
  for (size_t k = 0; k < READ_N_CATEGORIES; k++)
    out[n++] = r->written + k;
  for (size_t k = 0; k < READ_N_CATEGORIES; k++)
    out[n++] = r->filtered + k;
  assert(n == N_READ_COUNTS);
}

#define COUNTS_HEADER "flucut-counts 8"

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
  if (self->refs != NULL)
    for (size_t r = 0; r < self->refs->n_records; r++) {
      if (self->cover[r] == NULL)
        continue;
      size_t len = self->refs->records[r].seq_len;
      for (size_t side = 0; side < 2; side++)
        for (size_t cat = 0; cat < COVER_N_CATEGORIES; cat++) {
          /* one line per reference, side and category, with the nonzero
             position:count pairs */
          bool any = false;
          for (size_t pos = 0; pos <= len; pos++) {
            size_t v = self->cover[r][COVER_IDX(side, cat, pos, len)];
            if (v == 0)
              continue;
            if (!any) {
              fprintf(file, "cover %zu %zu %zu", r, side, cat);
              any = true;
            }
            fprintf(file, " %zu:%zu", pos, v);
          }
          if (any)
            fputc('\n', file);
        }
    }
  for (size_t g = 0; g < lengths_n_groups(self->lengths); g++)
    for (size_t b = 0; b < LENGTHS_N_BINS; b++) {
      size_t n = lengths_count(self->lengths, g, b);
      if (n > 0)
        fprintf(file, "length %zu %zu %zu\n", g, b, n);
    }
  cleavage_counts_save(&self->cleavage, file);
  marker_counts_save(&self->markers, file);
  subtype_reads_save(&self->subtype_reads, file);
  read_stats_save(self->stats, file);
  read_metadata_save(&self->metadata, file);
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
  if (sscanf(line, "cover %zu %zu %zu%n", &a, &b, &c, &n) == 3) {
    if (self->refs == NULL || a >= self->refs->n_records || b >= 2 ||
        c >= COVER_N_CATEGORIES)
      return false;
    size_t len = self->refs->records[a].seq_len;
    size_t *cov = ref_cover(self, a);
    rest = line + n;
    size_t pos, cnt;
    int m;
    while (sscanf(rest, " %zu:%zu%n", &pos, &cnt, &m) == 2) {
      if (pos > len)
        return false;
      cov[COVER_IDX(b, c, pos, len)] += cnt;
      rest += m;
    }
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
  if (strncmp(line, "subtype", 7) == 0)
    return subtype_reads_load_line(&self->subtype_reads, line);
  if (strncmp(line, "cleavage", 8) == 0)
    return cleavage_counts_load_line(&self->cleavage, line);
  if (strncmp(line, "marker-", 7) == 0)
    return marker_counts_load_line(&self->markers, line);
  if (line[0] == 'r')
    return read_stats_load_line(self->stats, line);
  if (line[0] == 'm')
    return read_metadata_load_line(&self->metadata, line);
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

/* add the interval [from, to) of the reference <ref>, clamped to it, to
   the coverage of the category <cat> */
static void add_span(counter_t *self, size_t ref, size_t cat, long from,
    long to) {
  long len = (long)self->refs->records[ref].seq_len;
  if (from < 0) from = 0;
  if (to > len) to = len;
  if (from >= to)
    return;
  size_t *c = ref_cover(self, ref);
  c[COVER_IDX(0, cat, (size_t)from, (size_t)len)]++;
  c[COVER_IDX(1, cat, (size_t)to, (size_t)len)]++;
}

/* add the span of a located read to the coverage, by its class (see
   read-cover.h); of a read with an internal deletion, only the two
   flanks */
static void count_cover(counter_t *self, const read_cover_t *cover) {
  if (!cover->located)
    return;
  if (cover->cover == COVER_DELETION) {
    add_span(self, cover->ref, COVER_DELETION, cover->from, cover->left_to);
    add_span(self, cover->ref, COVER_DELETION, cover->right_from,
             cover->to);
    return;
  }
  add_span(self, cover->ref, cover->cover, cover->from, cover->to);
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

void counter_add_located_read(counter_t *self, const amplicon_t *a,
    size_t len, const read_cover_t *cover) {
  assert(self != NULL);
  assert(a != NULL);
  assert(cover != NULL);
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
  if (cover->located)
    s->placed++;
  count_cover(self, cover);
  count_length(self, a);
}

void counter_add_read(counter_t *self, const amplicon_t *a, size_t len) {
  assert(self != NULL);
  read_cover_t cover;
  read_cover_locate(self->refs, NULL, a, &cover);
  counter_add_located_read(self, a, len, &cover);
}

void counter_add_split(counter_t *self, size_t n_parts, size_t len,
    size_t parts_len) {
  assert(self != NULL);
  assert(n_parts > 1);
  assert(parts_len <= len);
  counter_reads_t *s = &self->reads;
  assert(s->reads_in >= n_parts);
  s->split++;
  s->reads_in -= n_parts - 1;
  s->bases_in += len - parts_len;
}

void counter_add_input(counter_t *self, size_t len, const char *qual,
    const char *title, size_t title_len) {
  assert(self != NULL);
  read_stats_add(self->stats, COUNTER_STATS_INPUT, len, qual);
  if (title != NULL)
    read_metadata_add_title(&self->metadata, title, title_len);
}

void counter_add_cut(counter_t *self) {
  assert(self != NULL);
  self->reads.cut++;
}

void counter_add_cut_again(counter_t *self) {
  assert(self != NULL);
  self->reads.cut_again++;
}

void counter_add_written(counter_t *self, size_t len, size_t category,
    const char *qual) {
  assert(category < READ_N_CATEGORIES);
  self->reads.written[category]++;
  self->reads.reads_out++;
  self->reads.bases_out += len;
  read_stats_add(self->stats, COUNTER_STATS_WRITTEN, len, qual);
}

void counter_set_cleavage_reference(counter_t *self, const char *aa,
    const char *fusion) {
  assert(self != NULL);
  cleavage_counts_set_reference(&self->cleavage, aa, fusion);
}

void counter_add_cleavage_fusion(counter_t *self, const char *nt) {
  assert(self != NULL);
  cleavage_counts_add_fusion(&self->cleavage, nt);
}

void counter_add_cleavage_searched(counter_t *self) {
  assert(self != NULL);
  cleavage_counts_searched(&self->cleavage);
}

void counter_add_cleavage(counter_t *self, const char *nt) {
  assert(self != NULL);
  cleavage_counts_add(&self->cleavage, nt);
}

void counter_add_cleavage_control(counter_t *self, size_t i,
    const char *nt) {
  assert(self != NULL);
  cleavage_counts_add_control(&self->cleavage, i, nt);
}

void counter_set_subtype_reads(counter_t *self,
    const subtype_reads_t *subtype_reads) {
  assert(self != NULL);
  self->subtype_reads = *subtype_reads;
}

const subtype_reads_t* counter_subtype_reads(const counter_t *self) {
  assert(self != NULL);
  return &self->subtype_reads;
}

void counter_add_markers(counter_t *self,
    const marker_reference_t *const *chosen, long ref, const long *cuts,
    size_t n_cuts, const char *read, size_t len, double max_error_rate) {
  assert(self != NULL);
  marker_read_sites(chosen, ref, cuts, n_cuts, read, len, max_error_rate,
                    &self->markers);
}

const marker_counts_t* counter_markers(const counter_t *self) {
  assert(self != NULL);
  return &self->markers;
}

const cleavage_counts_t* counter_cleavage(const counter_t *self) {
  assert(self != NULL);
  return &self->cleavage;
}

const read_stats_t* counter_read_stats(const counter_t *self) {
  assert(self != NULL);
  return self->stats;
}

const read_metadata_t* counter_metadata(const counter_t *self) {
  assert(self != NULL);
  return &self->metadata;
}

void counter_add_too_short(counter_t *self) {
  self->reads.too_short++;
}

void counter_add_filtered(counter_t *self, size_t category) {
  assert(self != NULL);
  assert(category < READ_N_CATEGORIES);
  self->reads.filtered[category]++;
}

size_t counter_reads_filtered(const counter_reads_t *reads) {
  assert(reads != NULL);
  size_t n = 0;
  for (size_t k = 0; k < READ_N_CATEGORIES; k++)
    n += reads->filtered[k];
  return n;
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

bool counter_has_coverage(const counter_t *self, size_t ref) {
  assert(self != NULL);
  assert(self->refs != NULL && ref < self->refs->n_records);
  return self->cover[ref] != NULL;
}

size_t counter_coverage(const counter_t *self, size_t ref, size_t cat,
    size_t *out) {
  assert(self != NULL);
  assert(cat < COVER_N_CATEGORIES);
  assert(out != NULL);
  assert(self->refs != NULL && ref < self->refs->n_records);
  size_t len = self->refs->records[ref].seq_len;
  memset(out, 0, sizeof(size_t) * len);
  if (self->cover[ref] == NULL)
    return 0;
  const size_t *c = self->cover[ref];
  size_t running = 0, max = 0;
  for (size_t p = 0; p < len; p++) {
    running += c[COVER_IDX(0, cat, p, len)];
    running -= c[COVER_IDX(1, cat, p, len)];
    out[p] = running;
    if (running > max)
      max = running;
  }
  return max;
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

/* the coverage of reads with both ends mapped on a reference of 2000
   bases: one spanning it all, one with an internal deletion */
static void test_deletion_coverage(const scheme_t *sc) {
  static char seq[2001];
  memset(seq, 'A', 2000);
  fasta_record_t rec = {"r", "", seq, 2000};
  fasta_t refs = {&rec, 1};
  counter_t *c = counter_new(sc, &refs);
  amplicon_t a = test_amplicon(30, 2029, 0);
  a.left.origin = ORIGIN_TERMINAL;
  a.left.ref.ref_index = a.right.ref.ref_index = 0;
  a.right.ref.mapped = true;
  a.right.ref.ref_start = 1999;    /* the last base of the template */
  counter_add_read(c, &a, 2100);   /* 2000 bases over 2000: complete */
  a.keep_to = 429;
  counter_add_read(c, &a, 500);    /* 400 bases over 2000: a deletion */
  size_t *cov;
  XMALLOC(cov, sizeof(size_t) * 2000);
  test_eq_size(counter_coverage(c, 0, COVER_TERMINAL_PLUS, cov), 1);
  test_eq_size(cov[1000], 1);
  test_eq_size(cov[1999], 1);
  test_eq_size(counter_coverage(c, 0, COVER_DELETION, cov), 1);
  test_eq_size(cov[0], 1);
  test_eq_size(cov[199], 1);
  test_eq_size(cov[200], 0);
  test_eq_size(cov[1000], 0);
  test_eq_size(cov[1799], 0);
  test_eq_size(cov[1800], 1);
  test_eq_size(cov[1999], 1);
  /* with its junction located: each flank up to the junction */
  a.junction = (junction_t){.found = true, .left_to = 100,
                            .right_from = 1700};
  counter_add_read(c, &a, 500);
  test_eq_size(counter_coverage(c, 0, COVER_DELETION, cov), 2);
  test_eq_size(cov[99], 2);
  test_eq_size(cov[100], 1);
  test_eq_size(cov[1699], 0);
  test_eq_size(cov[1700], 1);
  test_eq_size(cov[1800], 2);
  free(cov);
  counter_delete(c);
}

/* the coverage of reads of which no end mapped, placed by a probe from
   inside them, on a reference of 2000 bases */
static void test_placed_coverage(const scheme_t *sc) {
  static char seq[2001];
  memset(seq, 'A', 2000);
  fasta_record_t rec = {"r", "", seq, 2000};
  fasta_t refs = {&rec, 1};
  counter_t *c = counter_new(sc, &refs);
  /* no primer: a fragment of unknown origin */
  amplicon_t a;
  memset(&a, 0, sizeof(a));
  a.strand = STRAND_NONE;
  a.keep_to = 299;
  a.placement = (placement_t){.found = true, .ref_index = 0,
                              .ref_len = 2000, .from = 500, .to = 800};
  counter_add_read(c, &a, 300);
  /* a 5' primer whose probe did not map: internal if the amplicon starts
     inside the reference, terminal at its end (here the end from which
     an amplicon on the minus strand of the reference starts) */
  a.strand = STRAND_PLUS;
  a.left.found = true;
  counter_add_read(c, &a, 300);
  a.placement.reverse = true;
  a.placement.from = 1690;
  a.placement.to = 1990;
  counter_add_read(c, &a, 300);
  size_t *cov;
  XMALLOC(cov, sizeof(size_t) * 2000);
  test_eq_size(counter_coverage(c, 0, COVER_NO_PRIMER, cov), 1);
  test_eq_size(cov[499], 0);
  test_eq_size(cov[500], 1);
  test_eq_size(cov[799], 1);
  test_eq_size(cov[800], 0);
  test_eq_size(counter_coverage(c, 0, COVER_INTERNAL_PLUS, cov), 1);
  test_eq_size(cov[500], 1);
  test_eq_size(counter_coverage(c, 0, COVER_TERMINAL_PLUS, cov), 1);
  test_eq_size(cov[1689], 0);
  test_eq_size(cov[1690], 1);
  test_eq_size(cov[1989], 1);
  free(cov);
  counter_delete(c);
}

/* the coverage of reads with one end mapped, on a reference of 2000
   bases: the span runs from the mapped end in the direction of the
   amplicon on the reference */
static void test_one_end_coverage(const scheme_t *sc) {
  static char seq[2001];
  memset(seq, 'A', 2000);
  fasta_record_t rec = {"r", "", seq, 2000};
  fasta_t refs = {&rec, 1};
  counter_t *c = counter_new(sc, &refs);
  size_t *cov;
  XMALLOC(cov, sizeof(size_t) * 2000);
  /* the left end, internal at 1000: forward, then on the minus strand */
  amplicon_t a = test_amplicon(30, 329, 1000);
  a.left.ref.ref_index = 0;
  a.left.ref.ref_len = 2000;
  counter_add_read(c, &a, 400);
  test_eq_size(counter_coverage(c, 0, COVER_INTERNAL_PLUS, cov), 1);
  test_eq_size(cov[999], 0);
  test_eq_size(cov[1000], 1);
  test_eq_size(cov[1299], 1);
  test_eq_size(cov[1300], 0);
  counter_clear(c);
  a.left.ref.reverse = true;
  counter_add_read(c, &a, 400);
  test_eq_size(counter_coverage(c, 0, COVER_INTERNAL_PLUS, cov), 1);
  test_eq_size(cov[700], 0);
  test_eq_size(cov[701], 1);
  test_eq_size(cov[1000], 1);
  test_eq_size(cov[1001], 0);
  counter_clear(c);
  /* the right end only, its probe forward at 100: an amplicon on the
     minus strand, which ends at 100 and started 300 bases further, far
     from the end of the reference (internal) */
  a.left.ref.mapped = false;
  a.right.ref.mapped = true;
  a.right.ref.ref_index = 0;
  a.right.ref.ref_len = 2000;
  a.right.ref.ref_start = 100;
  counter_add_read(c, &a, 400);
  test_eq_size(counter_coverage(c, 0, COVER_INTERNAL_PLUS, cov), 1);
  test_eq_size(cov[99], 0);
  test_eq_size(cov[100], 1);
  test_eq_size(cov[399], 1);
  test_eq_size(cov[400], 0);
  free(cov);
  counter_delete(c);
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
  counter_add_written(c, 100, READ_CATEGORY_INTERNAL, NULL);
  a = test_amplicon(30, 29, 7);
  counter_add_read(c, &a, 60);
  counter_add_too_short(c);
  counter_add_filtered(c, READ_CATEGORY_DELETION);
  const counter_reads_t *r = counter_reads(c);
  test_eq_size(r->reads_in, 2);
  test_eq_size(r->bases_in, 260);
  test_eq_size(r->reads_out, 1);
  test_eq_size(r->written[READ_CATEGORY_INTERNAL], 1);
  test_eq_size(r->written[READ_CATEGORY_FULL], 0);
  test_eq_size(r->n_ends[2], 2);
  test_eq_size(r->dimers, 1);
  test_eq_size(r->too_short, 1);
  test_eq_size(r->filtered[READ_CATEGORY_DELETION], 1);
  test_eq_size(counter_reads_filtered(r), 1);
  test_eq_size(r->origin[COUNTER_LEFT][ORIGIN_INTERNAL], 2);
  test_eq_size(r->origin[COUNTER_RIGHT][ORIGIN_TERMINAL], 2);
  /* both reads were placed on a reference by their ends */
  test_eq_size(r->placed, 2);
  test_eq_size(counter_primer(c, COUNTER_LEFT, 1), 2);
  test_eq_size(counter_primer(c, COUNTER_LEFT, 0), 0);
  test_eq_size(counter_primer(c, COUNTER_RIGHT, 0), 2);
  const size_t *s = counter_sites(c, 1, COUNTER_LEFT);
  test_assert(s != NULL);
  test_eq_size(s[7], 2);
  test_eq_size(s[6], 0);
  test_eq_size(counter_sites(c, 1, COUNTER_RIGHT)[7], 0);
  test_assert(counter_sites(c, 0, COUNTER_LEFT) == NULL);
  /* coverage: the read is internal on the plus strand of reference 1,
     its span [7, 7 + 100) clamped to the reference length (20) */
  size_t cov[32];
  test_eq_size(counter_coverage(c, 1, COVER_INTERNAL_PLUS, cov), 1);
  test_eq_size(cov[7], 1);
  test_eq_size(cov[19], 1);
  test_eq_size(cov[6], 0);
  test_eq_size(counter_coverage(c, 1, COVER_TERMINAL_PLUS, cov), 0);
  test_eq_size(counter_coverage(c, 0, COVER_INTERNAL_PLUS, cov), 0);
  /* the length of the template (100), in the group of reference 1; the
     dimer has none */
  test_eq_size(lengths_count(counter_lengths(c), 1, 10), 1);
  /* the sum of two counters */
  counter_merge(sum, c);
  counter_merge(sum, c);
  test_eq_size(counter_reads(sum)->reads_in, 4);
  test_eq_size(counter_reads(sum)->placed, 4);
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
  test_eq_size(counter_reads(loaded)->written[READ_CATEGORY_INTERNAL], 4);
  test_eq_size(counter_reads(loaded)->filtered[READ_CATEGORY_DELETION], 4);
  test_eq_size(counter_reads(loaded)->placed, 8);
  test_eq_size(counter_reads(loaded)->origin[COUNTER_LEFT][ORIGIN_INTERNAL],
               8);
  test_eq_size(counter_primer(loaded, COUNTER_LEFT, 1), 8);
  test_eq_size(counter_sites(loaded, 1, COUNTER_LEFT)[7], 8);
  test_eq_size(lengths_count(counter_lengths(loaded), 1, 10), 4);
  test_eq_size(counter_coverage(loaded, 1, COVER_INTERNAL_PLUS, cov), 4);
  test_eq_size(cov[7], 4);
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
  test_deletion_coverage(sc);
  test_placed_coverage(sc);
  test_one_end_coverage(sc);
  counter_clear(sum);
  test_eq_size(counter_reads(sum)->reads_in, 0);
  test_eq_size(counter_sites(sum, 1, COUNTER_LEFT)[7], 0);
  test_eq_size(counter_primer(sum, COUNTER_LEFT, 1), 0);
  counter_delete(c);
  counter_delete(sum);
  scheme_delete(sc);
  delete_logger(logger);
}
