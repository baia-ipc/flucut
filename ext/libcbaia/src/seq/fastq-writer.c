/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <string.h>
#include "fastq-writer.h"
#include "test.h"

void fastq_writer_record(strbuf_t *out, const char *title, size_t title_len,
    const char *seq, const char *qual, size_t len) {
  assert(out != NULL);
  assert(title != NULL);
  assert(len == 0 || (seq != NULL && qual != NULL));
  strbuf_append_char(out, '@');
  strbuf_append(out, title, title_len);
  strbuf_append_char(out, '\n');
  strbuf_append(out, seq, len);
  strbuf_append(out, "\n+\n", 3);
  strbuf_append(out, qual, len);
  strbuf_append_char(out, '\n');
}

/* --- unit tests --- */

void fastq_writer_test(void) {
  strbuf_t out;
  strbuf_init(&out);
  fastq_writer_record(&out, "r1 x=1", 6, "ACGT", "IIII", 4);
  fastq_writer_record(&out, "r2", 2, "", "", 0);
  const char *exp = "@r1 x=1\nACGT\n+\nIIII\n@r2\n\n+\n\n";
  test_eq(out.len, strlen(exp), size_t, "%zu");
  test_assert(memcmp(out.s, exp, out.len) == 0);
  strbuf_free(&out);
}
