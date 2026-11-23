/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <math.h>
#include <string.h>
#include "read-stats.h"
#include "test.h"

/* the error probability of each quality character (Phred+33) */
#define N_QUAL_CHARS 94
static double error_probability[N_QUAL_CHARS];
static bool probabilities_ready = false;

/* (set once, to the same values by any thread) */
static void init_probabilities(void) {
  for (int q = 0; q < N_QUAL_CHARS; q++)
    error_probability[q] = pow(10.0, -q / 10.0);
  probabilities_ready = true;
}

/* the quality of a read with the qualities <qual>: the Phred score of the
   mean error probability of its bases */
static double read_quality(const char *qual, size_t len) {
  if (!probabilities_ready)
    init_probabilities();
  double sum = 0;
  for (size_t i = 0; i < len; i++) {
    int q = qual[i] - '!';
    sum += error_probability[q < 0 ? 0 : q >= N_QUAL_CHARS ?
                             N_QUAL_CHARS - 1 : q];
  }
  /* (a little more: a read of equal qualities is counted in their bin,
     whatever the rounding of the sum) */
  return -10.0 * log10(sum / (double)len) + 1e-6;
}

void read_stats_clear(read_stats_t *self) {
  assert(self != NULL);
  memset(self, 0, sizeof(read_stats_t));
}

void read_stats_add(read_stats_t *self, size_t set, size_t len,
    const char *qual) {
  assert(self != NULL);
  assert(set < READ_STATS_N_SETS);
  read_stats_set_t *s = self->sets + set;
  if (s->reads == 0 || len < s->min_len)
    s->min_len = len;
  if (len > s->max_len)
    s->max_len = len;
  s->reads++;
  s->bases += len;
  size_t bin = len / READ_STATS_LEN_BIN;
  if (bin >= READ_STATS_LEN_BINS)
    bin = READ_STATS_LEN_BINS - 1;
  s->len_reads[bin]++;
  s->len_bases[bin] += len;
  if (qual == NULL || len == 0)
    return;
  double q = read_quality(qual, len);
  s->with_quality++;
  s->quality[q < 0 ? 0 : q >= READ_STATS_QUAL_MAX ? READ_STATS_QUAL_MAX :
             (size_t)q]++;
}

void read_stats_merge(read_stats_t *self, const read_stats_t *other) {
  assert(self != NULL);
  assert(other != NULL);
  for (size_t k = 0; k < READ_STATS_N_SETS; k++) {
    read_stats_set_t *s = self->sets + k;
    const read_stats_set_t *o = other->sets + k;
    if (o->reads == 0)
      continue;
    if (s->reads == 0 || o->min_len < s->min_len)
      s->min_len = o->min_len;
    if (o->max_len > s->max_len)
      s->max_len = o->max_len;
    s->reads += o->reads;
    s->bases += o->bases;
    s->with_quality += o->with_quality;
    for (size_t b = 0; b < READ_STATS_LEN_BINS; b++) {
      s->len_reads[b] += o->len_reads[b];
      s->len_bases[b] += o->len_bases[b];
    }
    for (size_t q = 0; q < READ_STATS_QUAL_BINS; q++)
      s->quality[q] += o->quality[q];
  }
}

/* the mean length of the reads of the bin of the lengths at which the
   running sum of <counts> (the reads or the bases of the bins, taken
   from the shortest reads if <ascending>, else from the longest) reaches
   half of <total> */
static size_t half_length(const read_stats_set_t *s, const size_t *counts,
    size_t total, bool ascending) {
  size_t sum = 0;
  for (size_t i = 0; i < READ_STATS_LEN_BINS; i++) {
    size_t b = ascending ? i : READ_STATS_LEN_BINS - 1 - i;
    sum += counts[b];
    if (counts[b] > 0 && 2 * sum >= total)
      return s->len_bases[b] / s->len_reads[b];
  }
  return 0;
}

void read_stats_summarize(const read_stats_t *self, size_t set,
    read_stats_summary_t *out) {
  assert(self != NULL);
  assert(set < READ_STATS_N_SETS);
  assert(out != NULL);
  const read_stats_set_t *s = self->sets + set;
  memset(out, 0, sizeof(read_stats_summary_t));
  out->reads = s->reads;
  out->bases = s->bases;
  if (s->reads == 0)
    return;
  out->min_len = s->min_len;
  out->max_len = s->max_len;
  out->mean_len = s->bases / s->reads;
  out->median_len = half_length(s, s->len_reads, s->reads, true);
  out->n50 = half_length(s, s->len_bases, s->bases, false);
  if (s->with_quality == 0)
    return;
  out->has_quality = true;
  size_t sum = 0, q10 = 0, q20 = 0;
  bool found = false;
  for (size_t q = 0; q < READ_STATS_QUAL_BINS; q++) {
    sum += s->quality[q];
    if (!found && 2 * sum >= s->with_quality) {
      /* (the middle of the bin) */
      out->median_quality = (double)q + 0.5;
      found = true;
    }
    if (q >= 10)
      q10 += s->quality[q];
    if (q >= 20)
      q20 += s->quality[q];
  }
  out->q10 = (double)q10 / (double)s->with_quality;
  out->q20 = (double)q20 / (double)s->with_quality;
}

const char* read_stats_set_name(size_t set) {
  assert(set < READ_STATS_N_SETS);
  return set == READ_STATS_INPUT ? "input" : "written";
}

void read_stats_save(const read_stats_t *self, FILE *file) {
  assert(self != NULL);
  assert(file != NULL);
  for (size_t k = 0; k < READ_STATS_N_SETS; k++) {
    const read_stats_set_t *s = self->sets + k;
    if (s->reads == 0)
      continue;
    fprintf(file, "rset %zu %zu %zu %zu\n", k, s->min_len, s->max_len,
            s->with_quality);
    for (size_t b = 0; b < READ_STATS_LEN_BINS; b++)
      if (s->len_reads[b] > 0)
        fprintf(file, "rlen %zu %zu %zu %zu\n", k, b, s->len_reads[b],
                s->len_bases[b]);
    for (size_t q = 0; q < READ_STATS_QUAL_BINS; q++)
      if (s->quality[q] > 0)
        fprintf(file, "rqual %zu %zu %zu\n", k, q, s->quality[q]);
  }
}

bool read_stats_load_line(read_stats_t *self, const char *line) {
  assert(self != NULL);
  assert(line != NULL);
  size_t k, a, b, c;
  if (sscanf(line, "rset %zu %zu %zu %zu", &k, &a, &b, &c) == 4) {
    if (k >= READ_STATS_N_SETS)
      return false;
    read_stats_set_t *s = self->sets + k;
    /* (the reads of the set are added by its lines of the lengths: until
       then, the lengths loaded are compared with those already there) */
    if (s->max_len == 0 || a < s->min_len)
      s->min_len = a;
    if (b > s->max_len)
      s->max_len = b;
    s->with_quality += c;
    return true;
  }
  if (sscanf(line, "rlen %zu %zu %zu %zu", &k, &a, &b, &c) == 4) {
    if (k >= READ_STATS_N_SETS || a >= READ_STATS_LEN_BINS)
      return false;
    read_stats_set_t *s = self->sets + k;
    s->len_reads[a] += b;
    s->len_bases[a] += c;
    s->reads += b;
    s->bases += c;
    return true;
  }
  if (sscanf(line, "rqual %zu %zu %zu", &k, &a, &b) == 3) {
    if (k >= READ_STATS_N_SETS || a >= READ_STATS_QUAL_BINS)
      return false;
    self->sets[k].quality[a] += b;
    return true;
  }
  return false;
}

void read_stats_test(void) {
  static read_stats_t stats, other;
  read_stats_summary_t s;
  read_stats_clear(&stats);
  read_stats_summarize(&stats, READ_STATS_INPUT, &s);
  test_eq_size(s.reads, 0);
  test_assert(!s.has_quality);
  /* three reads of 100, 300 and 1000 bases: half of the 1400 bases are
     in the longest one; the qualities: Q10 ('+'), Q20 ('5'), none */
  char qual[1001];
  memset(qual, '+', 100);
  read_stats_add(&stats, READ_STATS_INPUT, 100, qual);
  memset(qual, '5', 300);
  read_stats_add(&stats, READ_STATS_INPUT, 300, qual);
  read_stats_add(&stats, READ_STATS_INPUT, 1000, NULL);
  read_stats_summarize(&stats, READ_STATS_INPUT, &s);
  test_eq_size(s.reads, 3);
  test_eq_size(s.bases, 1400);
  test_eq_size(s.min_len, 100);
  test_eq_size(s.max_len, 1000);
  test_eq_size(s.mean_len, 466);
  test_eq_size(s.median_len, 300);
  test_eq_size(s.n50, 1000);
  test_assert(s.has_quality);
  test_assert(s.median_quality == 10.5);
  test_assert(s.q10 == 1.0 && s.q20 == 0.5);
  /* a bad base lowers the quality of a read more than the mean of the
     scores would: 9 bases of Q40 and one of Q1 give Q10 */
  read_stats_clear(&other);
  read_stats_add(&other, READ_STATS_WRITTEN, 10, "IIIIIIIII\"");
  read_stats_summarize(&other, READ_STATS_WRITTEN, &s);
  test_assert(s.median_quality == 10.5);
  /* merged, then saved and loaded */
  read_stats_merge(&other, &stats);
  read_stats_summarize(&other, READ_STATS_INPUT, &s);
  test_eq_size(s.reads, 3);
  test_eq_size(s.n50, 1000);
  char buf[4096], line[256];
  FILE *f = fmemopen(buf, sizeof(buf), "w+");
  read_stats_save(&other, f);
  rewind(f);
  read_stats_clear(&stats);
  while (fgets(line, sizeof(line), f) != NULL)
    test_assert(read_stats_load_line(&stats, line));
  fclose(f);
  test_assert(memcmp(&stats, &other, sizeof(read_stats_t)) == 0);
  test_assert(!read_stats_load_line(&stats, "rlen 9 0 1 1\n"));
}
