/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "logger.h"
#include "primers.h"
#include "scheme.h"
#include "strbuf.h"
#include "strutils.h"
#include "test.h"

/* a built-in scheme: the names of the primers of each end, NULL-ended */
#define SCHEME_MAX_PRIMERS 16
typedef struct {
  const char *name;
  const char *left[SCHEME_MAX_PRIMERS];
  const char *right[SCHEME_MAX_PRIMERS];
} builtin_scheme_t;

/* the built-in schemes, written by the build from
   metadata/primers.yaml (scripts/primers-table.awk) */
#define SCHEMES_TABLE_DEFINE
#include "primers-table.h"


static const char *builtin_names[SCHEMES_N_BUILTIN + 1];

const char* const* scheme_builtin_names(void) {
  for (size_t i = 0; i < SCHEMES_N_BUILTIN; i++)
    builtin_names[i] = SCHEMES_BUILTIN[i].name;
  builtin_names[SCHEMES_N_BUILTIN] = NULL;
  return builtin_names;
}

const char* scheme_default_for(virus_set_t set) {
  if (set == VIRUS_SET_A)
    return SCHEME_DEFAULT_A;
  if (set == VIRUS_SET_B)
    return SCHEME_DEFAULT_B;
  return SCHEME_DEFAULT_BOTH;
}

virus_set_t scheme_types(const scheme_t *self) {
  assert(self != NULL);
  virus_set_t set = 0;
  for (size_t i = 0; i < self->n_left; i++)
    set |= VIRUS_SET_OF(self->left[i].virus);
  for (size_t i = 0; i < self->n_right; i++)
    set |= VIRUS_SET_OF(self->right[i].virus);
  return set;
}

/* whether two names are equal, ignoring the case of ASCII letters */
static bool name_eq(const char *a, const char *b) {
  for (; *a != '\0' && *b != '\0'; a++, b++) {
    char ca = *a >= 'A' && *a <= 'Z' ? (char)(*a + 32) : *a;
    char cb = *b >= 'A' && *b <= 'Z' ? (char)(*b + 32) : *b;
    if (ca != cb)
      return false;
  }
  return *a == '\0' && *b == '\0';
}

/* the built-in scheme named <name> exactly (ignoring the case), NULL if
   there is none */
static const builtin_scheme_t* find_builtin(const char *name) {
  for (size_t i = 0; i < SCHEMES_N_BUILTIN; i++)
    if (name_eq(name, SCHEMES_BUILTIN[i].name))
      return SCHEMES_BUILTIN + i;
  return NULL;
}

/*
   The built-in schemes making up <name>: the scheme of that name, or
   the schemes named by its parts joined by '+' (e.g. opti+ibv-ga2), at
   most SCHEME_MAX_PARTS; their number, 0 if a part is not a built-in
   scheme
*/
#define SCHEME_MAX_PARTS 4
static size_t find_parts(const char *name,
    const builtin_scheme_t **parts) {
  const builtin_scheme_t *whole = find_builtin(name);
  if (whole != NULL) {
    parts[0] = whole;
    return 1;
  }
  size_t n = 0;
  const char *from = name;
  for (;;) {
    const char *plus = strchr(from, '+');
    size_t len = plus == NULL ? strlen(from) : (size_t)(plus - from);
    char part[64];
    if (n == SCHEME_MAX_PARTS || len == 0 || len >= sizeof(part))
      return 0;
    memcpy(part, from, len);
    part[len] = '\0';
    parts[n] = find_builtin(part);
    if (parts[n++] == NULL)
      return 0;
    if (plus == NULL)
      return n;
    from = plus + 1;
  }
}

bool scheme_is_builtin(const char *name) {
  assert(name != NULL);
  const builtin_scheme_t *parts[SCHEME_MAX_PARTS];
  return find_parts(name, parts) > 0;
}

/* the number of names of a NULL-ended list */
static size_t count_names(const char *const *names) {
  size_t n = 0;
  while (names[n] != NULL)
    n++;
  return n;
}

/* the primers of the lists of names of an end of the schemes <parts>
   (<right>: of the right end), appended to <dest>; their number (they
   are built in: a missing one is a programming error, not a user
   error) */
static size_t load_primers(const builtin_scheme_t **parts, size_t n_parts,
    bool right, primer_t *dest) {
  size_t n = 0;
  for (size_t k = 0; k < n_parts; k++) {
    const char *const *names = right ? parts[k]->right : parts[k]->left;
    for (size_t i = 0; names[i] != NULL; i++) {
      const primer_t *p = primers_builtin(names[i]);
      assert(p != NULL);
      dest[n++] = *p;
    }
  }
  return n;
}

/* the name of a scheme made of <parts>, as the names of the built-in
   schemes write it (newly allocated) */
static char* parts_name(const builtin_scheme_t **parts, size_t n_parts) {
  strbuf_t name;
  strbuf_init(&name);
  for (size_t k = 0; k < n_parts; k++) {
    if (k > 0)
      strbuf_append_char(&name, '+');
    strbuf_append_str(&name, parts[k]->name);
  }
  strbuf_append_char(&name, '\0');
  return name.s;
}

scheme_t* scheme_new_builtin(const char *name, Logger *logger) {
  assert(name != NULL);
  const builtin_scheme_t *parts[SCHEME_MAX_PARTS];
  size_t n_parts = find_parts(name, parts);
  if (n_parts == 0) {
    strbuf_t known;
    strbuf_init(&known);
    for (size_t i = 0; i < SCHEMES_N_BUILTIN; i++) {
      if (i > 0)
        strbuf_append_str(&known, ", ");
      strbuf_append_str(&known, SCHEMES_BUILTIN[i].name);
    }
    strbuf_append_char(&known, '\0');
    lprintf(logger, LOG_ERROR, "Unknown primer scheme \"%s\" (known "
            "schemes: %s, or several joined by +)\n", name, known.s);
    strbuf_free(&known);
    exit(EXIT_FAILURE);
  }
  size_t n_left = 0, n_right = 0;
  for (size_t k = 0; k < n_parts; k++) {
    n_left += count_names(parts[k]->left);
    n_right += count_names(parts[k]->right);
  }
  scheme_t *self;
  XMALLOC(self, sizeof(scheme_t));
  self->name = parts_name(parts, n_parts);
  XMALLOC(self->left, sizeof(primer_t) * (n_left > 0 ? n_left : 1));
  XMALLOC(self->right, sizeof(primer_t) * (n_right > 0 ? n_right : 1));
  self->n_left = load_primers(parts, n_parts, false, self->left);
  /* the primers of the right end, reverse-complemented: their sequences
     are stored one after the other in right_seqs */
  primer_t *forward;
  XMALLOC(forward, sizeof(primer_t) * (n_right > 0 ? n_right : 1));
  self->n_right = load_primers(parts, n_parts, true, forward);
  size_t total = 0;
  for (size_t i = 0; i < self->n_right; i++)
    total += forward[i].seq_len;
  XMALLOC(self->right_seqs, total > 0 ? total : 1);
  size_t offset = 0;
  for (size_t i = 0; i < self->n_right; i++) {
    primers_reverse(forward + i, self->right_seqs + offset,
                    self->right + i);
    offset += forward[i].seq_len;
  }
  free(forward);
  return self;
}

void scheme_replace(scheme_t *self, scheme_t *other) {
  assert(self != NULL);
  assert(other != NULL);
  scheme_t old = *self;
  *self = *other;
  *other = old;
  scheme_delete(other);
}

void scheme_delete(scheme_t *self) {
  if (self == NULL)
    return;
  free(self->name);
  free(self->left);
  free(self->right);
  free(self->right_seqs);
  free(self);
}

void scheme_test(void) {
  Logger *logger = new_logger(stderr, LOG_ERROR);
  scheme_t *s = scheme_new_builtin("mbt", logger);
  test_eq_str(s->name, "mbt");
  test_eq_size(s->n_left, 2);
  test_eq_size(s->n_right, 1);
  /* the primers of the left end, as they are ordered */
  test_eq_mem(s->left[0].seq, s->left[0].seq_len,
              "ACGCGTGATCAGCAAAAGCAGG");
  test_eq_mem(s->left[1].seq, s->left[1].seq_len,
              "ACGCGTGATCAGCRAAAGCAGG");
  test_eq_size(s->left[0].terminus_from, 10);
  /* the primer of the right end, reverse-complemented: the terminus is
     its prefix */
  test_eq_mem(s->right[0].seq, s->right[0].seq_len,
              "CCTTGTTTCTACTGATCACGCGT");
  test_eq_size(s->right[0].terminus_from, 0);
  test_eq_size(s->right[0].terminus_len, 13);
  scheme_delete(s);
  /* the scheme without tails */
  s = scheme_new_builtin("UNI", logger);
  test_eq_str(s->name, "uni");
  test_eq_mem(s->left[0].seq, s->left[0].seq_len, "AGCAAAAGCAGG");
  test_eq_mem(s->right[0].seq, s->right[0].seq_len, "CCTTGTTTCTACT");
  test_eq_size(s->right[0].terminus_len, 13);
  scheme_delete(s);
  /* the primers of the two ends can have different forms: an untailed
     Uni12 primer at the left end, a tailed one at the right end */
  /* the scheme of influenza B: the primers of the cocktail, each with
     the segments it binds; its right-end primers are stored
     reverse-complemented, with their terminus as their prefix */
  s = scheme_new_builtin(SCHEME_DEFAULT_B, logger);
  test_eq_size(s->n_left, 7);
  test_eq_size(s->n_right, 6);
  test_eq_int((int)scheme_types(s), (int)VIRUS_SET_B);
  test_eq_mem(s->left[0].seq, s->left[0].seq_len,
              "GGGGGGAGCAGAAGCGGAGC");
  test_eq_size(s->left[0].terminus_from, 6);
  test_eq_size(s->left[0].terminus_len, 14);
  test_eq_mem(s->right[0].seq, s->right[0].seq_len,
              "GCTCGTGTTTCTACTAATAACCCGG");
  test_eq_size(s->right[0].terminus_from, 0);
  test_eq_size(s->right[0].terminus_len, 15);
  test_assert(primers_binds_segment(s->right + 0, SEGMENT_PB1));
  test_assert(!primers_binds_segment(s->right + 0, SEGMENT_HA));
  scheme_delete(s);
  /* the default scheme holds the primers of both types */
  s = scheme_new_builtin(SCHEME_DEFAULT_BOTH, logger);
  test_eq_size(s->n_left, 8);
  test_eq_size(s->n_right, 7);
  test_eq_int((int)scheme_types(s), (int)VIRUS_SET_BOTH);
  scheme_delete(s);
  test_eq_str(scheme_default_for(VIRUS_SET_A), SCHEME_DEFAULT_A);
  test_eq_str(scheme_default_for(VIRUS_SET_B), SCHEME_DEFAULT_B);
  test_eq_str(scheme_default_for(VIRUS_SET_BOTH), SCHEME_DEFAULT_BOTH);
  s = scheme_new_builtin(SCHEME_DEFAULT_A, logger);
  test_eq_str(s->name, "uni12-mbt13");
  test_eq_size(s->n_left, 1);
  test_eq_mem(s->left[0].seq, s->left[0].seq_len, "AGCRAAAGCAGG");
  test_eq_size(s->left[0].terminus_from, 0);
  test_eq_size(s->left[0].terminus_len, 12);
  test_eq_mem(s->right[0].seq, s->right[0].seq_len,
              "CCTTGTTTCTACTGATCACGCGT");
  scheme_delete(s);
  /* a scheme made of two built-in ones: the primers of the first, then
     those of the second */
  test_assert(scheme_is_builtin("opti+ibv-ga2"));
  test_assert(scheme_is_builtin("UNI12-MBT13+IBV-GA2"));
  test_assert(!scheme_is_builtin("opti+"));
  test_assert(!scheme_is_builtin("opti+nothing"));
  test_assert(!scheme_is_builtin("auto"));
  s = scheme_new_builtin("Opti+ibv-ga2", logger);
  test_eq_str(s->name, "opti+ibv-ga2");
  test_eq_size(s->n_left, 9);
  test_eq_size(s->n_right, 7);
  test_eq_str(s->left[0].name, "opti1-f1");
  test_eq_str(s->left[2].name, "b-pbs-unif");
  test_eq_mem(s->right[0].seq, s->right[0].seq_len,
              "CCTTGTTTCTACTGGCGCGTAAC");
  test_eq_mem(s->right[1].seq, s->right[1].seq_len,
              "GCTCGTGTTTCTACTAATAACCCGG");
  test_eq_int((int)scheme_types(s), (int)VIRUS_SET_BOTH);
  /* replaced in place by another one */
  scheme_replace(s, scheme_new_builtin("uni", logger));
  test_eq_str(s->name, "uni");
  test_eq_size(s->n_left, 2);
  test_eq_size(s->n_right, 1);
  scheme_delete(s);
  /* the names of the built-in schemes: the default first, NULL-ended */
  const char *const *names = scheme_builtin_names();
  test_eq_str(names[0], SCHEME_DEFAULT);
  test_assert(names[SCHEMES_N_BUILTIN] == NULL);
  delete_logger(logger);
}
