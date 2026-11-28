/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <string.h>
#include <strings.h>
#include "test.h"
#include "virus.h"

static const char *type_names[VIRUS_N_TYPES] = {"A", "B"};

const char* virus_type_name(virus_type_t type) {
  assert(type == VIRUS_A || type == VIRUS_B);
  return type_names[type];
}

bool virus_type_parse(const char *s, virus_type_t *type) {
  assert(s != NULL);
  assert(type != NULL);
  /* the words which may precede the letter, so that "flu B",
     "influenza B" and "FluB" are accepted as well as "B" */
  static const char *prefixes[] = {"influenza", "flu", ""};
  for (size_t i = 0; i < sizeof(prefixes) / sizeof(prefixes[0]); i++) {
    size_t len = strlen(prefixes[i]);
    if (strncasecmp(s, prefixes[i], len) != 0)
      continue;
    const char *p = s + len;
    while (*p == ' ' || *p == '-' || *p == '_')
      p++;
    if (p[0] == '\0' || p[1] != '\0')
      continue;
    if (toupper((unsigned char)p[0]) == 'A') {
      *type = VIRUS_A;
      return true;
    }
    if (toupper((unsigned char)p[0]) == 'B') {
      *type = VIRUS_B;
      return true;
    }
  }
  return false;
}

/*
   Whether the letter <letter> names the type in <s> at <p>: as the word
   of a strain name ("B/Brisbane/60/2008") or of a seed name ("B_HA"),
   i.e. followed by '/' or '_' and not preceded by a letter or a digit
*/
static bool type_letter_at(const char *s, const char *p, char letter) {
  if (toupper((unsigned char)*p) != letter)
    return false;
  if (p > s && (isalnum((unsigned char)p[-1]) || p[-1] == '_'))
    return false;
  return p[1] == '/' || p[1] == '_';
}

bool virus_description_is_b(const char *description) {
  assert(description != NULL);
  /* "influenza B", "influenzavirus B", "influenza B virus" */
  for (const char *p = description; *p != '\0'; p++) {
    if (strncasecmp(p, "influenza", 9) != 0 ||
        (p > description && isalnum((unsigned char)p[-1])))
      continue;
    const char *q = p + 9;
    /* "influenzavirus B", "influenza virus B", "influenza B virus" */
    if (strncasecmp(q, "virus", 5) == 0)
      q += 5;
    while (*q == ' ' || *q == '-' || *q == '_')
      q++;
    if (strncasecmp(q, "virus", 5) == 0 &&
        !isalnum((unsigned char)q[5])) {
      q += 5;
      while (*q == ' ' || *q == '-' || *q == '_')
        q++;
    }
    if (toupper((unsigned char)*q) == 'B' && !isalnum((unsigned char)q[1]))
      return true;
    if (toupper((unsigned char)*q) == 'A' && !isalnum((unsigned char)q[1]))
      return false;
  }
  /* a strain name "B/..." or a seed name "B_..." */
  for (const char *p = description; *p != '\0'; p++) {
    if (type_letter_at(description, p, 'B'))
      return true;
    if (type_letter_at(description, p, 'A'))
      return false;
  }
  return false;
}

virus_type_t virus_type_of_description(const char *description) {
  return virus_description_is_b(description) ? VIRUS_B : VIRUS_A;
}

bool virus_set_has(virus_set_t set, virus_type_t type) {
  return (set & VIRUS_SET_OF(type)) != 0;
}

bool virus_set_parse(const char *s, virus_set_t *set) {
  assert(s != NULL);
  assert(set != NULL);
  if (strcasecmp(s, "both") == 0 || strcasecmp(s, "auto") == 0 ||
      strcasecmp(s, "ab") == 0 || strcasecmp(s, "a+b") == 0) {
    *set = VIRUS_SET_BOTH;
    return true;
  }
  virus_type_t type;
  if (!virus_type_parse(s, &type))
    return false;
  *set = VIRUS_SET_OF(type);
  return true;
}

const char* virus_set_name(virus_set_t set) {
  if (set == VIRUS_SET_A)
    return "A";
  if (set == VIRUS_SET_B)
    return "B";
  return "A and B";
}

virus_type_t virus_set_main(virus_set_t set) {
  return set == VIRUS_SET_B ? VIRUS_B : VIRUS_A;
}

void virus_test(void) {
  test_eq_str(virus_type_name(VIRUS_A), "A");
  test_eq_str(virus_type_name(VIRUS_B), "B");
  /* what the user may write */
  virus_type_t t = VIRUS_A;
  test_assert(virus_type_parse("B", &t));
  test_eq_int((int)t, (int)VIRUS_B);
  test_assert(virus_type_parse("b", &t));
  test_eq_int((int)t, (int)VIRUS_B);
  test_assert(virus_type_parse("influenza B", &t));
  test_eq_int((int)t, (int)VIRUS_B);
  test_assert(virus_type_parse("FluB", &t));
  test_eq_int((int)t, (int)VIRUS_B);
  test_assert(virus_type_parse("a", &t));
  test_eq_int((int)t, (int)VIRUS_A);
  test_assert(!virus_type_parse("auto", &t));
  test_assert(!virus_type_parse("C", &t));
  test_assert(!virus_type_parse("", &t));
  /* the type of a description */
  test_assert(virus_description_is_b("Influenza B virus "
      "(B/Brisbane/60/2008) segment 4 hemagglutinin (HA) gene"));
  test_assert(virus_description_is_b("B/Lee/1940 segment 1"));
  test_assert(virus_description_is_b("B_HA"));
  test_assert(virus_description_is_b("Influenzavirus B, complete genome"));
  test_assert(virus_description_is_b("Influenza B virus, segment 7"));
  test_assert(virus_description_is_b("influenza virus B segment 7"));
  test_assert(!virus_description_is_b("Influenza A virus "
      "(A/Puerto Rico/8/1934(H1N1)) segment 4 hemagglutinin (HA) gene"));
  test_assert(!virus_description_is_b("A_HA_H5"));
  test_assert(!virus_description_is_b("A/Hong Kong/1/68"));
  /* a description which names no type is taken for influenza A */
  test_assert(!virus_description_is_b("segment 4 hemagglutinin (HA)"));
  test_assert(!virus_description_is_b(""));
  /* the letter of a type is not a word of the strain name */
  test_assert(!virus_description_is_b("A/Bangkok/1/1979"));
  test_eq_int((int)virus_type_of_description("B_NS"), (int)VIRUS_B);
  test_eq_int((int)virus_type_of_description("A_NS"), (int)VIRUS_A);
  /* the sets of types */
  virus_set_t set = VIRUS_SET_A;
  test_assert(virus_set_parse("both", &set));
  test_eq_int((int)set, (int)VIRUS_SET_BOTH);
  test_assert(virus_set_has(set, VIRUS_A));
  test_assert(virus_set_has(set, VIRUS_B));
  test_assert(virus_set_parse("B", &set));
  test_eq_int((int)set, (int)VIRUS_SET_B);
  test_assert(!virus_set_has(set, VIRUS_A));
  test_assert(virus_set_has(set, VIRUS_B));
  test_assert(virus_set_parse("auto", &set));
  test_eq_int((int)set, (int)VIRUS_SET_BOTH);
  test_assert(!virus_set_parse("C", &set));
  test_eq_str(virus_set_name(VIRUS_SET_A), "A");
  test_eq_str(virus_set_name(VIRUS_SET_B), "B");
  test_eq_str(virus_set_name(VIRUS_SET_BOTH), "A and B");
  test_eq_int((int)virus_set_main(VIRUS_SET_B), (int)VIRUS_B);
  test_eq_int((int)virus_set_main(VIRUS_SET_BOTH), (int)VIRUS_A);
}
