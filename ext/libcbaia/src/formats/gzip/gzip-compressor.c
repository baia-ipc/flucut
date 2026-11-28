/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include "alloc.h"
#include "gzip-compressor.h"
#include "logger.h"
#include "strutils.h"
#include "test.h"

/* window bits for the gzip format */
#define GZIP_WINDOW_BITS (15 + 16)
#define GZIP_MEM_LEVEL 8

struct gzip_compressor_t {
  z_stream stream;
};

gzip_compressor_t* gzip_compressor_new(int level) {
  gzip_compressor_t *self;
  XCALLOC(self, 1, sizeof(gzip_compressor_t));
  if (deflateInit2(&self->stream, level, Z_DEFLATED, GZIP_WINDOW_BITS,
                   GZIP_MEM_LEVEL, Z_DEFAULT_STRATEGY) != Z_OK)
    fatal_printf("Initialization of the gzip compression failed\n");
  return self;
}

/* maximal amount of input or output of a deflate() call (avail_in and
   avail_out are unsigned int); can be lowered for testing */
#ifndef GZIP_CHUNK_MAX
#define GZIP_CHUNK_MAX UINT_MAX
#endif

static inline uInt chunk_size(size_t n) {
  return n > GZIP_CHUNK_MAX ? (uInt)GZIP_CHUNK_MAX : (uInt)n;
}

void gzip_compressor_run(gzip_compressor_t *self, const char *data,
    size_t len, strbuf_t *out) {
  assert(self != NULL);
  assert(data != NULL || len == 0);
  assert(out != NULL);
  z_stream *zs = &self->stream;
  size_t in_left = len;
  int ret;
  /* the input in chunks, so that it is not truncated by the conversion to
     uInt (4 GiB or more) */
  do {
    size_t room = deflateBound(zs, in_left);
    char *dest = strbuf_reserve(out, room);
    uInt in_chunk = chunk_size(in_left), out_chunk = chunk_size(room);
    zs->next_in = (Bytef *)(data + (len - in_left));
    zs->avail_in = in_chunk;
    zs->next_out = (Bytef *)dest;
    zs->avail_out = out_chunk;
    ret = deflate(zs, in_chunk == in_left ? Z_FINISH : Z_NO_FLUSH);
    if (ret != Z_OK && ret != Z_STREAM_END)
      fatal_printf("gzip compression failed\n");
    in_left -= in_chunk - zs->avail_in;
    strbuf_commit(out, out_chunk - zs->avail_out);
  } while (ret != Z_STREAM_END);
  deflateReset(zs);
}

void gzip_compressor_delete(gzip_compressor_t *self) {
  if (self == NULL)
    return;
  deflateEnd(&self->stream);
  free(self);
}

void gzip_compressor_test(const char *tmpdir) {
  gzip_compressor_t *gz = gzip_compressor_new(6);
  strbuf_t out;
  strbuf_init(&out);
  gzip_compressor_run(gz, "hello ", 6, &out);
  gzip_compressor_run(gz, "world\n", 6, &out);
  gzip_compressor_run(gz, "", 0, &out);
  char *path = strcat_alloc(tmpdir, "/cbaia-tests-gzip.gz");
  FILE *f = fopen(path, "wb");
  test_assert(f != NULL);
  fwrite(out.s, 1, out.len, f);
  fclose(f);
  gzFile in = gzopen(path, "rb");
  test_assert(in != NULL);
  char buf[64];
  int n = gzread(in, buf, sizeof(buf));
  gzclose(in);
  test_eq(n, 12, int, "%d");
  test_assert(memcmp(buf, "hello world\n", 12) == 0);
  /* a larger block (several chunks if GZIP_CHUNK_MAX is lowered) */
  size_t big_len = 100000;
  char *big, *back;
  XMALLOC(big, big_len);
  XMALLOC(back, big_len + 1);
  for (size_t i = 0; i < big_len; i++)
    big[i] = "ACGT\n"[(i * 7 + i / 13) % 5];
  strbuf_clear(&out);
  gzip_compressor_run(gz, big, big_len, &out);
  f = fopen(path, "wb");
  test_assert(f != NULL);
  fwrite(out.s, 1, out.len, f);
  fclose(f);
  in = gzopen(path, "rb");
  test_assert(in != NULL);
  n = gzread(in, back, (unsigned)big_len + 1);
  gzclose(in);
  test_eq(n, (int)big_len, int, "%d");
  test_assert(memcmp(big, back, big_len) == 0);
  free(big);
  free(back);
  unlink(path);
  free(path);
  strbuf_free(&out);
  gzip_compressor_delete(gz);
}
