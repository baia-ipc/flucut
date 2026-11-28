/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdlib.h>
#include "alloc.h"
#include "json-coverage.h"
#include "strutils.h"

/* the names of the classes of reads (COVER_*) */
static const char *const names[COVER_N_CATEGORIES] = {
  "terminal_plus", "terminal_minus", "internal_plus", "internal_minus",
  "no_primer", "not_located", "deletion"};

bool json_coverage_has(const counter_t *c, size_t ref,
    const fasta_t *refs) {
  if (!counter_has_coverage(c, ref))
    return false;
  /* (reads can be counted for a reference without any depth on it: a
     span outside the reference) */
  size_t *depth, max = 0;
  XMALLOC(depth, sizeof(size_t) * (refs->records[ref].seq_len + 1));
  for (size_t k = 0; k < COVER_N_CATEGORIES && max == 0; k++)
    max = counter_coverage(c, ref, k, depth);
  free(depth);
  return max > 0;
}

void json_coverage_write(FILE *f, const fasta_t *refs, const counter_t *c,
    size_t ref) {
  assert(f != NULL);
  assert(refs != NULL);
  assert(ref < refs->n_records);
  size_t len = refs->records[ref].seq_len;
  size_t n = len < JSON_COVERAGE_POINTS ? len : JSON_COVERAGE_POINTS;
  fputs("{\"reference\": ", f);
  fputs_json_string(refs->records[ref].name, f);
  fprintf(f, ", \"length\": %zu, \"positions\": [", len);
  for (size_t i = 0; i < n; i++)
    fprintf(f, "%s%zu", i > 0 ? ", " : "",
            n > 1 ? i * (len - 1) / (n - 1) : 0);
  fputs("], \"depth\": {", f);
  size_t *depth;
  XMALLOC(depth, sizeof(size_t) * (len > 0 ? len : 1));
  bool first = true;
  for (size_t k = 0; k < COVER_N_CATEGORIES; k++) {
    if (counter_coverage(c, ref, k, depth) == 0)
      continue;
    fprintf(f, "%s\"%s\": [", first ? "" : ", ", names[k]);
    first = false;
    for (size_t i = 0; i < n; i++)
      fprintf(f, "%s%zu", i > 0 ? ", " : "",
              depth[n > 1 ? i * (len - 1) / (n - 1) : 0]);
    fputc(']', f);
  }
  free(depth);
  fputs("}}", f);
}
