/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include "json-reads.h"
#include "strutils.h"

void json_reads_lengths(FILE *f, const counter_t *c) {
  assert(f != NULL);
  fputc('{', f);
  for (size_t k = 0; k < counter_read_stats(c)->n_sets; k++) {
    read_stats_summary_t s;
    read_stats_summarize(counter_read_stats(c), k, &s);
    const char *name = read_stats_set_name(counter_read_stats(c), k);
    fprintf(f, "%s\"%s\": {\"reads\": %zu, \"bases\": %zu, \"min\": %zu, "
            "\"max\": %zu, \"mean\": %zu, \"median\": %zu, \"n50\": %zu}",
            k > 0 ? ", " : "", name, s.reads, s.bases,
            s.min_len, s.max_len, s.mean_len, s.median_len, s.n50);
  }
  fputc('}', f);
}

void json_reads_quality(FILE *f, const counter_t *c) {
  assert(f != NULL);
  const read_stats_t *stats = counter_read_stats(c);
  fputc('{', f);
  for (size_t k = 0; k < stats->n_sets; k++) {
    read_stats_summary_t s;
    read_stats_summarize(stats, k, &s);
    fprintf(f, "%s\"%s\": ", k > 0 ? ", " : "", read_stats_set_name(stats, k));
    if (!s.has_quality)
      fputs("null", f);
    else
      fprintf(f, "{\"reads\": %zu, \"median\": %.1f, \"q10\": %.4f, "
              "\"q20\": %.4f}", stats->sets[k].with_quality,
              s.median_quality, s.q10, s.q20);
  }
  fputc('}', f);
}

void json_reads_sequencing(FILE *f, const counter_t *c) {
  assert(f != NULL);
  const read_metadata_t *m = counter_metadata(c);
  fputc('{', f);
  for (size_t i = 0; i < READ_METADATA_N; i++) {
    const char *text = read_metadata_text(m, i);
    fprintf(f, "%s\"%s\": ", i > 0 ? ", " : "", read_metadata_key(i));
    if (text == NULL)
      fputs("null", f);
    else
      fputs_json_string(text, f);
  }
  fputc('}', f);
}
