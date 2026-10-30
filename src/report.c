/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "logger.h"
#include "report.h"
#include "segments.h"
#include "test.h"

const char* report_end_name(size_t end) {
  assert(end < COUNTER_N_ENDS);
  return end == COUNTER_LEFT ? "left" : "right";
}

const char* report_origin_name(end_origin_t origin) {
  static const char *names[] = {"unknown", "terminal", "internal"};
  return names[origin];
}

/* order of the sites: decreasing reads, then reference, end, position */
static int site_cmp(const void *a, const void *b) {
  const report_site_t *x = a, *y = b;
  if (x->reads != y->reads)
    return x->reads > y->reads ? -1 : 1;
  if (x->ref != y->ref)
    return x->ref < y->ref ? -1 : 1;
  if (x->end != y->end)
    return x->end < y->end ? -1 : 1;
  return x->pos < y->pos ? -1 : x->pos > y->pos ? 1 : 0;
}

/* add the sites of the counts <s> of the positions of a reference of
   length <len> (end <end>) to <sites> */
static report_site_t* add_sites(report_site_t *sites, size_t *n,
    size_t *cap, const size_t *s, size_t len, size_t ref, size_t end) {
  report_site_t site = {0};
  bool open = false;
  for (size_t p = 0; p < len; p++) {
    if (s[p] == 0)
      continue;
    if (open && p - (site.to - 1) > REPORT_SITE_GAP) {
      XCHECKSPACEINC(sites, cap, sizeof(report_site_t) * (*n + 1), 64);
      sites[(*n)++] = site;
      open = false;
    }
    if (!open) {
      site = (report_site_t){.ref = ref, .end = end, .from = p + 1,
                             .pos = p + 1};
      open = true;
    } else if (s[p] > s[site.pos - 1])
      site.pos = p + 1;
    site.to = p + 1;
    site.reads += s[p];
  }
  if (open) {
    XCHECKSPACEINC(sites, cap, sizeof(report_site_t) * (*n + 1), 64);
    sites[(*n)++] = site;
  }
  return sites;
}

report_site_t* report_sites(const fasta_t *refs, const counter_t *c,
    size_t *n) {
  assert(c != NULL);
  assert(n != NULL);
  report_site_t *sites = NULL;
  size_t cap = 0;
  *n = 0;
  if (refs == NULL)
    return NULL;
  for (size_t r = 0; r < refs->n_records; r++)
    for (size_t e = 0; e < COUNTER_N_ENDS; e++) {
      const size_t *s = counter_sites(c, r, e);
      if (s != NULL)
        sites = add_sites(sites, n, &cap, s, refs->records[r].seq_len, r,
                          e);
    }
  if (*n > 1)
    qsort(sites, *n, sizeof(report_site_t), site_cmp);
  return sites;
}

void report_lengths(const fasta_t *refs, const counter_t *c, size_t group,
    report_lengths_t *out) {
  assert(c != NULL);
  assert(out != NULL);
  const lengths_t *l = counter_lengths(c);
  assert(group < lengths_n_groups(l));
  out->ref_len = refs != NULL && group < refs->n_records ?
                 refs->records[group].seq_len : 0;
  out->name = refs == NULL ? "all" : group < refs->n_records ?
              refs->records[group].name : "not_located";
  if (refs != NULL && group < refs->n_records)
    segments_label(out->name, refs->records[group].description, out->label,
                   sizeof(out->label));
  else
    snprintf(out->label, sizeof(out->label), "%s", out->name);
  /* (no short amplicons without a reference: expected length 0) */
  lengths_summarize(l, group, out->ref_len, &out->summary);
}

size_t report_short_complete(const fasta_t *refs, const counter_t *c) {
  if (refs == NULL)
    return 0;
  size_t n = 0;
  for (size_t g = 0; g < refs->n_records; g++) {
    report_lengths_t l;
    report_lengths(refs, c, g, &l);
    n += l.summary.short_count;
  }
  return n;
}

/* --- unit tests --- */

void report_test(void) {
  Logger *logger = new_logger(stderr, LOG_ERROR);
  scheme_t *sc = scheme_new_builtin(SCHEME_DEFAULT, logger);
  char seq[1001];
  memset(seq, 'A', 1000);
  seq[1000] = '\0';
  fasta_record_t recs[2] = {{"r0", "", seq, 1000}, {"r1", "", seq, 1000}};
  fasta_t refs = {recs, 2};
  counter_t *c = counter_new(sc, &refs);
  size_t n;
  test_assert(report_sites(&refs, c, &n) == NULL);
  test_eq_size(n, 0);
  /* left ends of reference 1 at 0-based 100 (3 reads), 103 (5), 107 (1):
     one site, 101-108, most at 104; at 300 (2): another; a right end of
     reference 0 at 100 (4) */
  amplicon_t a;
  memset(&a, 0, sizeof(a));
  a.left.found = true;
  a.left.origin = ORIGIN_INTERNAL;
  a.left.ref.mapped = true;
  a.left.ref.ref_index = 1;
  long left[] = {100, 100, 100, 103, 103, 103, 103, 103, 107, 300, 300};
  for (size_t i = 0; i < sizeof(left) / sizeof(left[0]); i++) {
    a.left.ref.ref_start = left[i];
    counter_add_read(c, &a, 50);
  }
  a.left.found = false;
  a.right = a.left;
  a.right.found = true;
  a.right.ref.ref_index = 0;
  a.right.ref.ref_start = 100;
  for (size_t i = 0; i < 4; i++)
    counter_add_read(c, &a, 50);
  report_site_t *s = report_sites(&refs, c, &n);
  test_eq_size(n, 3);
  test_eq_size(s[0].ref, 1);
  test_eq_size(s[0].end, COUNTER_LEFT);
  test_eq_size(s[0].reads, 9);
  test_eq_size(s[0].pos, 104);
  test_eq_size(s[0].from, 101);
  test_eq_size(s[0].to, 108);
  test_eq_size(s[1].ref, 0);
  test_eq_size(s[1].end, COUNTER_RIGHT);
  test_eq_size(s[1].reads, 4);
  test_eq_size(s[1].pos, 101);
  test_eq_size(s[2].reads, 2);
  test_eq_size(s[2].pos, 301);
  free(s);
  /* the complete amplicons: none */
  report_lengths_t l;
  report_lengths(&refs, c, 2, &l);
  test_eq_str(l.name, "not_located");
  test_eq_size(l.ref_len, 0);
  test_eq_size(l.summary.total, 0);
  report_lengths(&refs, c, 1, &l);
  test_eq_str(l.name, "r1");
  test_eq_size(l.ref_len, 1000);
  test_eq_size(report_short_complete(&refs, c), 0);
  test_eq_size(report_short_complete(NULL, c), 0);
  test_eq_str(report_end_name(COUNTER_RIGHT), "right");
  test_eq_str(report_origin_name(ORIGIN_INTERNAL), "internal");
  counter_delete(c);
  scheme_delete(sc);
  delete_logger(logger);
}
