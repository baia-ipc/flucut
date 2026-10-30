/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "adapters.h"
#include "adaptor-search.h"
#include "alloc.h"
#include "test.h"

/* a barcode of unknown sequence (see adapters.h) */
#define BC "NNNNNNNNNNNNNNNNNNNNNNNN"

/* the kits (see adapters.h): the adapter at the start of a read and that
   at its end (as they are found in the read; NULL: none) */
static const struct {
  const char *name;
  const char *start, *end;
} kits[] = {
  {"nb", "ATTGCTAAGGTTAA" BC "CAGCACCT", "GGTGCTG" BC "TTAACCTTAGCAAT"},
  {"rbk", "GCTTGGGTGTTTAACC" BC
          "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA", NULL},
  {"lsk", "CCTGTACTTCGTTCAGTTACGTATTGC", "AGCAATACGTAACTGAAC"},
};
#define N_KITS (sizeof(kits) / sizeof(kits[0]))

struct adapters_t {
  bool used[N_KITS];
  long window;
  adaptor_search_t *search;
};

/* the kit named <name> (length <len>), N_KITS if none */
static size_t kit_of_name(const char *name, size_t len) {
  for (size_t i = 0; i < N_KITS; i++)
    if (strlen(kits[i].name) == len && strncmp(kits[i].name, name, len) == 0)
      return i;
  return N_KITS;
}

/* parse <list> into <used>; false if a name is unknown */
static bool parse_kits(const char *list, bool used[N_KITS]) {
  for (size_t i = 0; i < N_KITS; i++)
    used[i] = false;
  if (strcmp(list, "none") == 0)
    return true;
  const char *p = list;
  while (true) {
    const char *comma = strchr(p, ',');
    size_t len = comma == NULL ? strlen(p) : (size_t)(comma - p);
    size_t k = kit_of_name(p, len);
    if (k == N_KITS)
      return false;
    used[k] = true;
    if (comma == NULL)
      return true;
    p = comma + 1;
  }
}

bool adapters_valid(const char *kits_list) {
  bool used[N_KITS];
  return parse_kits(kits_list, used);
}

adapters_t* adapters_new(const char *kits_list, long window,
    Logger *logger) {
  adapters_t *self;
  XMALLOC(self, sizeof(adapters_t));
  LCHECK(parse_kits(kits_list, self->used), logger,
         "Unknown adapter kit in '%s' (kits: " ADAPTERS_KITS ", or none)\n",
         kits_list);
  self->window = window;
  self->search = adaptor_search_new();
  return self;
}

bool adapters_any(const adapters_t *self) {
  for (size_t i = 0; i < N_KITS; i++)
    if (self->used[i])
      return true;
  return false;
}

const char* adapters_kit_name(size_t kit) {
  assert(kit < N_KITS);
  return kits[kit].name;
}

void adapters_delete(adapters_t *self) {
  if (self == NULL)
    return;
  adaptor_search_delete(self->search);
  free(self);
}

/* the errors allowed in <adapter> (see adapters.h) */
static long max_errors(const char *adapter) {
  size_t fixed = 0;
  for (const char *p = adapter; *p != '\0'; p++)
    if (*p != 'N')
      fixed++;
  return (long)(ADAPTERS_ERROR_RATE * (double)fixed);
}

/* search the adapters at one end of <seq> (see adapters.h) */
static void search(adapters_t *self, const char *seq, size_t len,
    bool at_start, adapter_hit_t *hit) {
  hit->found = false;
  for (size_t k = 0; k < N_KITS; k++) {
    const char *adapter = at_start ? kits[k].start : kits[k].end;
    if (!self->used[k] || adapter == NULL)
      continue;
    size_t adapter_len = strlen(adapter);
    /* the region of the read in which the adapter may lie */
    size_t region_len = (size_t)self->window + adapter_len;
    if (region_len > len)
      region_len = len;
    size_t offset = at_start ? 0 : len - region_len;
    /* an adapter with a barcode may not overhang the read: the part
       left in the read could be mostly the barcode, which matches any
       sequence */
    bool partial = strchr(adapter, 'N') == NULL;
    adaptor_hit_t h;
    if (!adaptor_search_run(self->search, adapter, adapter_len,
                            seq + offset, region_len, max_errors(adapter),
                            partial && at_start, partial && !at_start, &h))
      continue;
    long start = h.start + (long)offset, end = h.end + (long)offset;
    if (at_start ? start > self->window :
        (long)len - 1 - end > self->window)
      continue;
    /* the innermost (see adapters.h) */
    if (hit->found && (at_start ? end <= hit->end : start >= hit->start))
      continue;
    hit->found = true;
    hit->kit = k;
    hit->start = start;
    hit->end = end;
    hit->n_errors = h.n_errors;
  }
}

void adapters_search_start(adapters_t *self, const char *seq, size_t len,
    adapter_hit_t *hit) {
  assert(self != NULL);
  assert(hit != NULL);
  search(self, seq, len, true, hit);
}

void adapters_search_end(adapters_t *self, const char *seq, size_t len,
    adapter_hit_t *hit) {
  assert(self != NULL);
  assert(hit != NULL);
  search(self, seq, len, false, hit);
}

void adapters_test(void) {
  test_assert(adapters_valid("nb"));
  test_assert(adapters_valid("nb,rbk,lsk"));
  test_assert(adapters_valid("none"));
  test_assert(!adapters_valid("nb,xx"));
  test_assert(!adapters_valid(""));
  adapters_t *a = adapters_new("nb,lsk", 100, NULL);
  test_assert(adapters_any(a));
  adapter_hit_t hit;
  /* flank, ligation adapter, native barcode NB01, template */
  const char *read = "GGCATTACCTGTACTTCGTTCAGTTACGTATTGCTAAGGTTAA"
                     "CACAAAGACACCGACAACTTTCTTCAGCACCTGGAATTACTTAACTAGC";
  adapters_search_start(a, read, strlen(read), &hit);
  test_assert(hit.found);
  test_eq_str(adapters_kit_name(hit.kit), "nb");
  test_eq_long(hit.end, 74L);
  /* without the barcode kit: only the ligation adapter */
  adapters_t *l = adapters_new("lsk", 100, NULL);
  adapters_search_start(l, read, strlen(read), &hit);
  test_assert(hit.found);
  test_eq_long(hit.end, 33L);
  adapters_delete(l);
  /* template, the reverse complement of the bottom barcode adapter of
     NB01 (with a substitution), flank */
  const char *end = "GGAATTACTTAACTAGCGGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTT"
                    "AACCTTAGGAATTTACC";
  adapters_search_end(a, end, strlen(end), &hit);
  test_assert(hit.found);
  test_eq_long(hit.start, 17L);
  test_eq_long(hit.n_errors, 1L);
  /* no adapter */
  const char *none = "GGAATTACTTAACTAGCAAAATGGAAACAATATCACTAATAACTATACTA";
  adapters_search_start(a, none, strlen(none), &hit);
  test_assert(!hit.found);
  adapters_search_end(a, none, strlen(none), &hit);
  test_assert(!hit.found);
  /* outside the window */
  adapters_t *w = adapters_new("nb", 10, NULL);
  adapters_search_start(w, read, strlen(read), &hit);
  test_assert(!hit.found);
  adapters_delete(w);
  adapters_t *n = adapters_new("none", 100, NULL);
  test_assert(!adapters_any(n));
  adapters_search_start(n, read, strlen(read), &hit);
  test_assert(!hit.found);
  adapters_delete(n);
  adapters_delete(a);
}
