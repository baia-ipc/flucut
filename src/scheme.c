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
typedef struct {
  const char *name;
  const char *left[4];
  const char *right[4];
} builtin_scheme_t;

static const builtin_scheme_t builtin[] = {
  {"uni12-mbt13", {"uni12r", NULL}, {"mbtuni-13", NULL}},
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
  s = scheme_new_builtin(SCHEME_DEFAULT, logger);
  test_eq_str(s->name, "uni12-mbt13");
  test_eq_size(s->n_left, 1);
  test_eq_mem(s->left[0].seq, s->left[0].seq_len, "AGCRAAAGCAGG");
  test_eq_size(s->left[0].terminus_from, 0);
  test_eq_size(s->left[0].terminus_len, 12);
  test_eq_mem(s->right[0].seq, s->right[0].seq_len,
              "CCTTGTTTCTACTGATCACGCGT");
  scheme_delete(s);
  const char *const *names = scheme_builtin_names();
  test_eq_str(names[0], "uni12-mbt13");
  test_assert(names[N_BUILTIN] == NULL);
  delete_logger(logger);
}
