/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include "segments.h"
#include "test.h"

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
    return n;
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

void segments_label(const char *name, const char *description, char *buf,
    size_t size) {
  int s = description == NULL ? 0 : segments_of_description(description);
  if (s == 0)
    snprintf(buf, size, "%s", name);
  else
    snprintf(buf, size, "segment %d (%s), %s", s, segments_gene(s), name);
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
}
