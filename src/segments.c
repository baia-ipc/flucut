/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include "segments.h"
#include "test.h"
#include "virus.h"

static const char *genes[SEGMENTS_N] = {"PB2", "PB1", "PA", "HA", "NP",
                                        "NA", "M", "NS"};

/* the words which name the gene of a segment (see segments.h) */
static const struct {
  const char *word;
  int segment;
} gene_words[] = {
  {"PB2", 1}, {"PB1", 2}, {"PA", 3}, {"HA", 4}, {"hemagglutinin", 4},
  {"NP", 5}, {"nucleoprotein", 5}, {"nucleocapsid", 5}, {"NA", 6},
  {"neuraminidase", 6}, {"M", 7}, {"M1", 7}, {"M2", 7}, {"MP", 7},
  {"matrix", 7}, {"NS", 8}, {"NS1", 8}, {"NS2", 8}, {"NEP", 8},
  {"nonstructural", 8},
  /* "non-structural": the next word is checked */
  {"non", -8},
};
#define N_GENE_WORDS (sizeof(gene_words) / sizeof(gene_words[0]))

/*
   The segment of the number <number> in <type> and the number of a
   segment: influenza B numbers PB1 1 and PB2 2, the other way round
   from influenza A, and numbers every other segment alike (see
   segments.h). The exchange is its own inverse.
*/
static int exchange_pb(int n) {
  return n == 1 ? 2 : n == 2 ? 1 : n;
}

int segments_number(int segment, virus_type_t type) {
  assert(segment >= 1 && segment <= SEGMENTS_N);
  return type == VIRUS_B ? exchange_pb(segment) : segment;
}

int segments_of_number(int number, virus_type_t type) {
  assert(number >= 1 && number <= SEGMENTS_N);
  return type == VIRUS_B ? exchange_pb(number) : number;
}

/* "segment <n>" in <s>, 0 if none */
static int segment_number(const char *s) {
  static const char key[] = "segment";
  const size_t key_len = sizeof(key) - 1;
  for (const char *p = s; *p != '\0'; p++) {
    if (strncasecmp(p, key, key_len) != 0 ||
        (p > s && isalnum((unsigned char)p[-1])))
      continue;
    const char *q = p + key_len;
    while (*q == ' ' || *q == ':' || *q == '_')
      q++;
    if (*q >= '1' && *q <= '0' + SEGMENTS_N &&
        !isalnum((unsigned char)q[1]))
      return *q - '0';
  }
  return 0;
}

/* the segment of the word <w> (length <len>; <next>: the text after it),
   0 if none */
static int segment_of_word(const char *w, size_t len, const char *next) {
  for (size_t i = 0; i < N_GENE_WORDS; i++) {
    const char *g = gene_words[i].word;
    /* the gene names (upper case) exactly, the protein names in any
       case */
    bool any_case = islower((unsigned char)g[0]);
    if (strlen(g) != len || (any_case ? strncasecmp(w, g, len) :
                             strncmp(w, g, len)) != 0)
      continue;
    if (gene_words[i].segment > 0)
      return gene_words[i].segment;
    /* "non-structural", "non structural" */
    if ((*next == '-' || *next == ' ') &&
        strncasecmp(next + 1, "structural", 10) == 0)
      return -gene_words[i].segment;
  }
  return 0;
}

int segments_of_description(const char *description) {
  assert(description != NULL);
  int n = segment_number(description);
  if (n > 0)
    return segments_of_number(n,
        virus_type_of_description(description));
  const char *p = description;
  while (*p != '\0') {
    while (*p != '\0' && !isalnum((unsigned char)*p))
      p++;
    const char *w = p;
    while (isalnum((unsigned char)*p))
      p++;
    int s = p > w ? segment_of_word(w, (size_t)(p - w), p) : 0;
    if (s > 0)
      return s;
  }
  return 0;
}

const char *segments_gene(int segment) {
  assert(segment >= 1 && segment <= SEGMENTS_N);
  return genes[segment - 1];
}

int segments_parse(const char *s, virus_type_t type) {
  assert(s != NULL);
  if (s[0] >= '1' && s[0] <= '0' + SEGMENTS_N && s[1] == '\0')
    return segments_of_number(s[0] - '0', type);
  for (size_t i = 0; i < N_GENE_WORDS; i++)
    if (gene_words[i].segment > 0 && strcasecmp(s, gene_words[i].word) == 0)
      return gene_words[i].segment;
  return 0;
}

void segments_name(int segment, virus_type_t type, char *buf,
    size_t size) {
  snprintf(buf, size, "segment %d (%s)", segments_number(segment, type),
           segments_gene(segment));
}

void segments_name_in(int segment, const char *description, char *buf,
    size_t size) {
  segments_name(segment, description == NULL ? VIRUS_A :
                virus_type_of_description(description), buf, size);
}

void segments_label(const char *name, const char *description, char *buf,
    size_t size) {
  int s = description == NULL ? 0 : segments_of_description(description);
  if (s == 0)
    snprintf(buf, size, "%s", name);
  else
    snprintf(buf, size, "segment %d (%s), %s",
             segments_number(s, virus_type_of_description(description)),
             segments_gene(s), name);
}

void segments_test(void) {
  test_eq_int(segments_of_description("Influenza A virus (A/Puerto Rico/"
          "8/1934(H1N1)) segment 4 hemagglutinin (HA) gene"), 4);
  test_eq_int(segments_of_description("Influenza A virus ha gene for "
          "Hemagglutinin, genomic RNA, strain A/Hong Kong/1073/99(H9N2)"),
          4);
  test_eq_int(segments_of_description("A_HA_H5"), 4);
  test_eq_int(segments_of_description("A_NA_N1"), 6);
  test_eq_int(segments_of_description("A_MP"), 7);
  test_eq_int(segments_of_description("matrix protein 2 (M2) and matrix "
          "protein 1 (M1) genes"), 7);
  test_eq_int(segments_of_description("non-structural protein"), 8);
  test_eq_int(segments_of_description("polymerase PB1 (PB1) gene"), 2);
  test_eq_int(segments_of_description("SEGMENT:6"), 6);
  /* a numbered segment of influenza B is read in its own numbering */
  test_eq_int(segments_of_description("Influenza B virus "
          "(B/Lee/1940) segment 1, complete sequence"), 2);
  test_eq_int(segments_of_description("Influenza B virus "
          "(B/Lee/1940) segment 2, complete sequence"), 1);
  test_eq_int(segments_of_description("Influenza B virus "
          "(B/Lee/1940) segment 4, complete sequence"), 4);
  test_eq_int(segments_of_description("B_PB1"), 2);
  test_eq_int(segments_of_description("B_HA"), 4);
  test_eq_int(segments_of_description("B_MP"), 7);
  test_eq_int(segments_of_description("segment 9 something"), 0);
  test_eq_int(segments_of_description("strain A/Hong Kong/1/68"), 0);
  test_eq_int(segments_of_description("na gene"), 0);
  test_eq_int(segments_of_description(""), 0);
  char buf[64];
  segments_label("NC_004908.1", "ha gene for Hemagglutinin", buf,
                 sizeof(buf));
  test_eq_str(buf, "segment 4 (HA), NC_004908.1");
  segments_label("ref1", "", buf, sizeof(buf));
  test_eq_str(buf, "ref1");
  segments_label("NC_002204.1", "Influenza B virus (B/Lee/1940) "
                 "segment 1, complete sequence", buf, sizeof(buf));
  test_eq_str(buf, "segment 1 (PB1), NC_002204.1");
  test_eq_int(segments_parse("4", VIRUS_A), 4);
  test_eq_int(segments_parse("HA", VIRUS_A), 4);
  test_eq_int(segments_parse("ha", VIRUS_A), 4);
  test_eq_int(segments_parse("Neuraminidase", VIRUS_A), 6);
  test_eq_int(segments_parse("MP", VIRUS_A), 7);
  test_eq_int(segments_parse("9", VIRUS_A), 0);
  test_eq_int(segments_parse("0", VIRUS_A), 0);
  test_eq_int(segments_parse("44", VIRUS_A), 0);
  test_eq_int(segments_parse("non", VIRUS_A), 0);
  test_eq_int(segments_parse("", VIRUS_A), 0);
  /* the numbers of influenza B: its segment 1 is PB1, its 2 is PB2;
     the gene names are the same in both types */
  test_eq_int(segments_parse("1", VIRUS_B), 2);
  test_eq_int(segments_parse("2", VIRUS_B), 1);
  test_eq_int(segments_parse("3", VIRUS_B), 3);
  test_eq_int(segments_parse("PB1", VIRUS_B), 2);
  test_eq_int(segments_number(1, VIRUS_A), 1);
  test_eq_int(segments_number(1, VIRUS_B), 2);
  test_eq_int(segments_number(2, VIRUS_B), 1);
  test_eq_int(segments_number(4, VIRUS_B), 4);
  test_eq_int(segments_of_number(segments_number(2, VIRUS_B), VIRUS_B), 2);
  segments_name(4, VIRUS_A, buf, sizeof(buf));
  test_eq_str(buf, "segment 4 (HA)");
  segments_name(1, VIRUS_B, buf, sizeof(buf));
  test_eq_str(buf, "segment 2 (PB1)");
  segments_name_in(2, "Influenza B virus (B/Lee/1940) segment 1", buf,
                   sizeof(buf));
  test_eq_str(buf, "segment 1 (PB1)");
  segments_name_in(2, NULL, buf, sizeof(buf));
  test_eq_str(buf, "segment 2 (PB1)");
}
