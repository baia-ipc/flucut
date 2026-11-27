/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "reference.h"
#include "segments.h"
#include "subtype-reads.h"
#include "test.h"

/* the size of the name of a subtype in a description (e.g. "H10N11") */
#define SUBTYPE_NAME_SIZE 16

void subtype_reads_clear(subtype_reads_t *self) {
  assert(self != NULL);
  memset(self, 0, sizeof(*self));
}

const char* subtype_reads_gene(size_t which) {
  return which == SUBTYPE_READS_HA ? "HA" : "NA";
}

char subtype_reads_letter(size_t which) {
  return which == SUBTYPE_READS_HA ? 'H' : 'N';
}

/* the segment (segments.h) of <which> */
static int segment_of(size_t which) {
  return which == SUBTYPE_READS_HA ? SEGMENT_HA : SEGMENT_NA;
}

int subtype_reads_number(const char *description, int segment) {
  char name[SUBTYPE_NAME_SIZE];
  if (segment != SEGMENT_HA && segment != SEGMENT_NA)
    return 0;
  if (!reference_description_subtype(description, name, sizeof(name)))
    return 0;
  /* (the lineages of influenza B are not subtypes) */
  const char *n = strchr(name, 'N');
  if (name[0] != 'H' || n == NULL)
    return 0;
  int number = atoi(segment == SEGMENT_HA ? name + 1 : n + 1);
  return number >= 1 && number <= SUBTYPE_READS_MAX_NUMBER ? number : 0;
}

int* subtype_reads_numbers(const fasta_t *fa, const int *segment) {
  int *number;
  XMALLOC(number, sizeof(int) * (fa->n_records > 0 ? fa->n_records : 1));
  for (size_t i = 0; i < fa->n_records; i++)
    number[i] = subtype_reads_number(fa->records[i].description,
                                     segment[i]);
  return number;
}

/* add the probes of a read which vote in the segment <s> to <counts> */
static void add_segment(subtype_reads_segment_t *counts, int s,
    const int *segment, const int *number, const size_t *read_votes,
    size_t n) {
  size_t k = 0;
  int first = 0;
  bool agree = true;
  for (size_t j = 0; j < n; j++) {
    size_t ref = read_votes[j];
    if (segment[ref] != s || number[ref] == 0)
      continue;
    counts->votes[number[ref]]++;
    if (k == 0)
      first = number[ref];
    else if (number[ref] != first)
      agree = false;
    k++;
  }
  if (k < 2)
    return;
  counts->reads[k]++;
  if (agree)
    counts->agreeing[first]++;
}

void subtype_reads_add(subtype_reads_t *self, const int *segment,
    const int *number, const size_t *read_votes, size_t n) {
  assert(self != NULL);
  if (n > SUBTYPE_READS_MAX_PROBES)
    n = SUBTYPE_READS_MAX_PROBES;
  for (size_t which = 0; which < SUBTYPE_READS_N_SEGMENTS; which++)
    add_segment(self->segments + which, segment_of(which), segment, number,
                read_votes, n);
}

void subtype_reads_merge(subtype_reads_t *self,
    const subtype_reads_t *other) {
  assert(self != NULL);
  assert(other != NULL);
  self->samples += other->samples;
  for (size_t which = 0; which < SUBTYPE_READS_N_SEGMENTS; which++) {
    subtype_reads_segment_t *a = self->segments + which;
    const subtype_reads_segment_t *b = other->segments + which;
    for (int g = 0; g <= SUBTYPE_READS_MAX_NUMBER; g++) {
      a->votes[g] += b->votes[g];
      a->agreeing[g] += b->agreeing[g];
    }
    for (size_t k = 0; k <= SUBTYPE_READS_MAX_PROBES; k++)
      a->reads[k] += b->reads[k];
  }
}

/* the reads expected to be all for a subtype with the share <share> of
   the votes, if the probes voted independently of their read */
static double expected_reads(const subtype_reads_segment_t *counts,
    double share) {
  double expected = 0, power = share;
  for (size_t k = 2; k <= SUBTYPE_READS_MAX_PROBES; k++) {
    power *= share;
    expected += (double)counts->reads[k] * power;
  }
  return expected;
}

bool subtype_reads_above_threshold(double excess, size_t reads) {
  return excess >= SUBTYPE_READS_MIN_EXCESS &&
         excess >= SUBTYPE_READS_MIN_SHARE * (double)reads;
}

/* the second subtype of <counts>, whose leading one is <out->leading>
   and whose votes are <total>: the other subtype with reads whose
   excess is largest */
static void second_subtype(const subtype_reads_segment_t *counts,
    size_t total, subtype_reads_evidence_t *out) {
  for (int g = 1; g <= SUBTYPE_READS_MAX_NUMBER; g++) {
    if (g == out->leading || counts->agreeing[g] == 0)
      continue;
    double expected = expected_reads(counts,
        (double)counts->votes[g] / (double)total);
    double excess = (double)counts->agreeing[g] - expected;
    if (out->second != 0 && excess <= out->second_excess)
      continue;
    out->second = g;
    out->second_reads = counts->agreeing[g];
    out->second_expected = expected;
    out->second_excess = excess;
  }
}

void subtype_reads_evidence(const subtype_reads_t *self, size_t which,
    subtype_reads_evidence_t *out) {
  assert(self != NULL);
  assert(which < SUBTYPE_READS_N_SEGMENTS);
  const subtype_reads_segment_t *counts = self->segments + which;
  memset(out, 0, sizeof(*out));
  size_t total = 0;
  for (int g = 1; g <= SUBTYPE_READS_MAX_NUMBER; g++) {
    total += counts->votes[g];
    if (counts->votes[g] > counts->votes[out->leading])
      out->leading = g;
  }
  for (size_t k = 2; k <= SUBTYPE_READS_MAX_PROBES; k++)
    out->reads += counts->reads[k];
  out->known = self->samples == 1 && out->reads > 0 && out->leading != 0;
  if (!out->known)
    return;
  out->leading_reads = counts->agreeing[out->leading];
  second_subtype(counts, total, out);
  if (out->second == 0)
    return;
  out->second_share = out->second_excess / (double)out->reads;
  out->above_threshold =
      subtype_reads_above_threshold(out->second_excess, out->reads);
}

void subtype_reads_log(const subtype_reads_t *self, const char *file,
    Logger *logger) {
  for (size_t which = 0; which < SUBTYPE_READS_N_SEGMENTS; which++) {
    subtype_reads_evidence_t e;
    subtype_reads_evidence(self, which, &e);
    if (!e.known)
      continue;
    char letter = subtype_reads_letter(which);
    const char *gene = subtype_reads_gene(which);
    if (e.second == 0) {
      lprintf(logger, LOG_INFO, "%s subtypes in the reads of %s: %zu "
              "reads with several voting probes, %zu all for %c%d, none "
              "for another subtype\n", gene, file, e.reads,
              e.leading_reads, letter, e.leading);
      continue;
    }
    lprintf(logger, LOG_INFO, "%s subtypes in the reads of %s: %zu reads "
            "with several voting probes, %zu all for %c%d, %zu all for "
            "%c%d (%.1f expected if the probes voted independently of "
            "their read; excess %.1f, %.2f%% of the reads%s)\n", gene,
            file, e.reads, e.leading_reads, letter, e.leading,
            e.second_reads, letter, e.second, e.second_expected,
            e.second_excess, 100 * e.second_share,
            e.above_threshold ? "; above the threshold" : "");
  }
}

void subtype_reads_save(const subtype_reads_t *self, FILE *file) {
  assert(self != NULL);
  assert(file != NULL);
  if (self->samples == 0)
    return;
  fprintf(file, "subtype-samples %zu\n", self->samples);
  for (size_t which = 0; which < SUBTYPE_READS_N_SEGMENTS; which++) {
    const subtype_reads_segment_t *counts = self->segments + which;
    for (int g = 1; g <= SUBTYPE_READS_MAX_NUMBER; g++)
      if (counts->votes[g] > 0)
        fprintf(file, "subtype-votes %zu %d %zu %zu\n", which, g,
                counts->votes[g], counts->agreeing[g]);
    for (size_t k = 2; k <= SUBTYPE_READS_MAX_PROBES; k++)
      if (counts->reads[k] > 0)
        fprintf(file, "subtype-reads %zu %zu %zu\n", which, k,
                counts->reads[k]);
  }
}

bool subtype_reads_load_line(subtype_reads_t *self, const char *line) {
  assert(self != NULL);
  size_t which, a, b, c;
  if (sscanf(line, "subtype-samples %zu", &a) == 1) {
    self->samples += a;
    return true;
  }
  if (sscanf(line, "subtype-votes %zu %zu %zu %zu", &which, &a, &b, &c)
      == 4) {
    if (which >= SUBTYPE_READS_N_SEGMENTS || a < 1 ||
        a > SUBTYPE_READS_MAX_NUMBER)
      return false;
    self->segments[which].votes[a] += b;
    self->segments[which].agreeing[a] += c;
    return true;
  }
  if (sscanf(line, "subtype-reads %zu %zu %zu", &which, &a, &b) == 3) {
    if (which >= SUBTYPE_READS_N_SEGMENTS || a < 2 ||
        a > SUBTYPE_READS_MAX_PROBES)
      return false;
    self->segments[which].reads[a] += b;
    return true;
  }
  return false;
}

/* the numbers of the subtypes, read from the descriptions */
static void test_numbers(void) {
  const char *h5n1 = "Influenza A virus (A/goose/X/1/2023(H5N1)) segment 4";
  test_eq(subtype_reads_number(h5n1, SEGMENT_HA), 5, int, "%d");
  test_eq(subtype_reads_number(h5n1, SEGMENT_NA), 1, int, "%d");
  test_eq(subtype_reads_number(h5n1, 1), 0, int, "%d");
  test_eq(subtype_reads_number("A/duck/X/2/1999(H10N11)", SEGMENT_NA), 11,
          int, "%d");
  test_eq(subtype_reads_number("A/duck/X/2/1999(H10N11)", SEGMENT_HA), 10,
          int, "%d");
  test_eq(subtype_reads_number("Influenza A virus segment 4", SEGMENT_HA),
          0, int, "%d");
  test_eq(subtype_reads_number("Influenza B virus (B/X/1/2008) segment 4 "
                               "Victoria lineage", SEGMENT_HA), 0, int,
          "%d");
  test_eq(subtype_reads_number(NULL, SEGMENT_HA), 0, int, "%d");
}

/* add <n> reads with the votes <read_votes> */
static void add_n(subtype_reads_t *self, const int *segment,
    const int *number, const size_t *read_votes, size_t n_votes, size_t n) {
  for (size_t i = 0; i < n; i++)
    subtype_reads_add(self, segment, number, read_votes, n_votes);
}

/* the references of the tests: two of H5 (0, 1), one of H3 (2), one of
   N1 (3), one of N2 (4), one of another segment (5) */
static const int test_segment[] = {4, 4, 4, 6, 6, 1};
static const int test_number[] = {5, 5, 3, 1, 2, 0};

/* a sample of one virus: the probes of its reads go to the two
   references of its subtype, a few single probes to another one */
static void test_one_virus(void) {
  subtype_reads_t c;
  subtype_reads_clear(&c);
  c.samples = 1;
  const size_t h5[] = {0, 1, 0, 3, 3, 5}, stray[] = {0, 2, 1, 0};
  const size_t single[] = {0};
  add_n(&c, test_segment, test_number, h5, 6, 200);
  add_n(&c, test_segment, test_number, stray, 4, 10);
  add_n(&c, test_segment, test_number, single, 1, 5);
  subtype_reads_evidence_t e;
  subtype_reads_evidence(&c, SUBTYPE_READS_HA, &e);
  test_assert(e.known);
  test_eq(e.leading, 5, int, "%d");
  test_eq(e.second, 0, int, "%d");
  test_eq_size(e.reads, 210);
  test_eq_size(e.leading_reads, 200);
  test_eq_size(e.second_reads, 0);
  test_assert(!e.above_threshold);
  test_eq_size(c.segments[SUBTYPE_READS_HA].votes[5], 635);
  test_eq_size(c.segments[SUBTYPE_READS_HA].votes[3], 10);
  /* the NA: the reads with two probes in it, all N1 */
  subtype_reads_evidence(&c, SUBTYPE_READS_NA, &e);
  test_eq(e.leading, 1, int, "%d");
  test_eq_size(e.reads, 200);
  test_eq_size(e.leading_reads, 200);
}

/* a sample with the reads of a second HA subtype */
static void test_two_subtypes(void) {
  subtype_reads_t c;
  subtype_reads_clear(&c);
  c.samples = 1;
  const size_t h5[] = {0, 1, 0}, h3[] = {2, 2};
  add_n(&c, test_segment, test_number, h5, 3, 190);
  add_n(&c, test_segment, test_number, h3, 2, 10);
  subtype_reads_evidence_t e;
  subtype_reads_evidence(&c, SUBTYPE_READS_HA, &e);
  test_eq(e.leading, 5, int, "%d");
  test_eq(e.second, 3, int, "%d");
  test_eq_size(e.reads, 200);
  test_eq_size(e.second_reads, 10);
  /* the share of H3 is 20 of 590 votes: 190 reads with 3 probes and 10
     with 2 are expected to be all H3 in 0.019 cases */
  test_assert(e.second_expected > 0.018 && e.second_expected < 0.020);
  test_assert(e.second_excess > 9.98 && e.second_excess < 9.99);
  test_assert(e.second_share > 0.0499 && e.second_share < 0.05);
  test_assert(e.above_threshold);
  /* nothing in the NA: no read with several probes there */
  subtype_reads_evidence(&c, SUBTYPE_READS_NA, &e);
  test_assert(!e.known);
  /* the counts of two samples together say nothing */
  subtype_reads_t sum;
  subtype_reads_clear(&sum);
  subtype_reads_merge(&sum, &c);
  subtype_reads_merge(&sum, &c);
  test_eq_size(sum.samples, 2);
  test_eq_size(sum.segments[SUBTYPE_READS_HA].agreeing[3], 20);
  subtype_reads_evidence(&sum, SUBTYPE_READS_HA, &e);
  test_assert(!e.known);
  test_assert(!e.above_threshold);
}

/* the threshold: enough reads, and a large enough share of the reads */
static void test_threshold(void) {
  test_assert(!subtype_reads_above_threshold(
      SUBTYPE_READS_MIN_EXCESS - 0.1, 100));
  test_assert(subtype_reads_above_threshold(SUBTYPE_READS_MIN_EXCESS,
                                            100));
  size_t many = (size_t)(2 * SUBTYPE_READS_MIN_EXCESS /
                         SUBTYPE_READS_MIN_SHARE);
  test_assert(!subtype_reads_above_threshold(SUBTYPE_READS_MIN_EXCESS,
                                             many));
  test_assert(subtype_reads_above_threshold(
      2 * SUBTYPE_READS_MIN_EXCESS, many));
}

/* the counts saved and loaded are the same */
static void test_save(void) {
  subtype_reads_t c, loaded;
  subtype_reads_clear(&c);
  subtype_reads_clear(&loaded);
  c.samples = 1;
  const size_t h5[] = {0, 1, 0, 3, 4}, h3[] = {2, 2, 4, 4};
  add_n(&c, test_segment, test_number, h5, 5, 7);
  add_n(&c, test_segment, test_number, h3, 4, 2);
  char *buf = NULL;
  size_t size = 0;
  FILE *f = open_memstream(&buf, &size);
  test_assert(f != NULL);
  subtype_reads_save(&c, f);
  fclose(f);
  for (char *line = strtok(buf, "\n"); line != NULL;
       line = strtok(NULL, "\n"))
    test_assert(subtype_reads_load_line(&loaded, line));
  free(buf);
  test_assert(memcmp(&c, &loaded, sizeof(c)) == 0);
  test_assert(!subtype_reads_load_line(&loaded, "subtype-votes 2 5 1 1"));
  test_assert(!subtype_reads_load_line(&loaded, "subtype-reads 0 1 4"));
  test_assert(!subtype_reads_load_line(&loaded, "reads 1 2 3"));
}

void subtype_reads_test(void) {
  test_numbers();
  test_one_virus();
  test_two_subtypes();
  test_threshold();
  test_save();
}
