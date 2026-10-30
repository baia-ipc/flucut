/* (c) 2026 Giorgio Gonnella */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bgzf-compressor.h"
#include "file.h"
#include "test.h"
#include "test-bam.h"
#include "ubam-writer.h"

/* compression level of the test data */
#define TEST_BAM_COMPRESSLEVEL 6

void test_bam_data(strbuf_t *out, const char *text,
    const char *const *titles, size_t n_titles) {
  strbuf_t data;
  strbuf_init(&data);
  if (text != NULL)
    ubam_writer_header(&data, text);
  for (size_t i = 0; i < n_titles; i++)
    ubam_writer_record(&data, titles[i], strlen(titles[i]), "ACGT", "IIII",
                       4);
  bgzf_compressor_t *c = bgzf_compressor_new(TEST_BAM_COMPRESSLEVEL);
  bgzf_compressor_run(c, data.s, data.len, out);
  bgzf_compressor_delete(c);
  strbuf_free(&data);
}

char *test_bam_file(const char *tmpdir, const char *name, const char *text,
    const char *const *titles, size_t n_titles, bool eof) {
  strbuf_t out;
  strbuf_init(&out);
  test_bam_data(&out, text, titles, n_titles);
  if (eof)
    strbuf_append(&out, (const char *)bgzf_eof, BGZF_EOF_LEN);
  char *path = path_join(tmpdir, name);
  FILE *f = fopen(path, "wb");
  test_assert(f != NULL);
  test_eq_size(fwrite(out.s, 1, out.len, f), out.len);
  fclose(f);
  strbuf_free(&out);
  return path;
}
