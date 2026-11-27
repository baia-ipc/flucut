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

/* the primers of influenza B of each end (the IBV-GA2 cocktail, see
   primers.h), written once as they are used by two schemes */
#define IBV_GA2_LEFT "b-pbs-unif", "b-pa-unif", "b-hana-unif", \
  "b-np-unif", "b-m-uni3f", "b-mg-uni3f", "b-ns-uni3f"
#define IBV_GA2_RIGHT "b-pbs-unir", "b-pa-unir", "b-hana-unir", \
  "b-np-unir", "b-m-uni3r", "b-ns-uni3r"

static const builtin_scheme_t builtin[] = {
  {SCHEME_DEFAULT_BOTH, {"uni12r", IBV_GA2_LEFT, NULL},
   {"mbtuni-13", IBV_GA2_RIGHT, NULL}},
  {SCHEME_DEFAULT_A, {"uni12r", NULL}, {"mbtuni-13", NULL}},
  {SCHEME_DEFAULT_B, {IBV_GA2_LEFT, NULL}, {IBV_GA2_RIGHT, NULL}},
  {"mbt", {"mbtuni-12", "mbtuni-12.4r", NULL}, {"mbtuni-13", NULL}},
  {"uni", {"uni12", "uni12r", NULL}, {"uni13", NULL}},
};
#define N_BUILTIN ((size_t)(sizeof(builtin) / sizeof(builtin[0])))

static const char *builtin_names[N_BUILTIN + 1];

const char* const* scheme_builtin_names(void) {
  for (size_t i = 0; i < N_BUILTIN; i++)
    builtin_names[i] = builtin[i].name;
  builtin_names[N_BUILTIN] = NULL;
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

/* the number of names of a NULL-ended list */
static size_t count_names(const char *const *names) {
  size_t n = 0;
  while (names[n] != NULL)
    n++;
  return n;
}

/* the primers of a list of names (they are built in: a missing one is a
   programming error, not a user error) */
static void load_primers(const char *const *names, size_t n,
    primer_t *dest) {
  for (size_t i = 0; i < n; i++) {
    const primer_t *p = primers_builtin(names[i]);
    assert(p != NULL);
    dest[i] = *p;
  }
}

scheme_t* scheme_new_builtin(const char *name, Logger *logger) {
  assert(name != NULL);
  const builtin_scheme_t *b = NULL;
  for (size_t i = 0; i < N_BUILTIN && b == NULL; i++)
    if (name_eq(name, builtin[i].name))
      b = builtin + i;
  if (b == NULL) {
    strbuf_t known;
    strbuf_init(&known);
    for (size_t i = 0; i < N_BUILTIN; i++) {
      if (i > 0)
        strbuf_append_str(&known, ", ");
      strbuf_append_str(&known, builtin[i].name);
    }
    strbuf_append_char(&known, '\0');
    lprintf(logger, LOG_ERROR, "Unknown primer scheme \"%s\" (known "
            "schemes: %s)\n", name, known.s);
    strbuf_free(&known);
    exit(EXIT_FAILURE);
  }
  scheme_t *self;
  XMALLOC(self, sizeof(scheme_t));
  self->name = xstrdup(b->name);
  self->n_left = count_names(b->left);
  self->n_right = count_names(b->right);
  XMALLOC(self->left, sizeof(primer_t) * self->n_left);
  XMALLOC(self->right, sizeof(primer_t) * self->n_right);
  load_primers(b->left, self->n_left, self->left);
  /* the primers of the right end, reverse-complemented: their sequences
     are stored one after the other in right_seqs */
  primer_t *forward;
  XMALLOC(forward, sizeof(primer_t) * self->n_right);
  load_primers(b->right, self->n_right, forward);
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
  /* the names of the built-in schemes: the default first, NULL-ended */
  const char *const *names = scheme_builtin_names();
  test_eq_str(names[0], SCHEME_DEFAULT);
  test_assert(names[N_BUILTIN] == NULL);
  delete_logger(logger);
}
