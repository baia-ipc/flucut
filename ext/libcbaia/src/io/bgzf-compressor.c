/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include "alloc.h"
#include "bgzf-compressor.h"
#include "logger.h"
#include "strutils.h"
#include "test.h"

/* maximum uncompressed size of a block (as in htslib), so that the
   compressed block is not larger than 64 KiB */
#define BGZF_BLOCK_DATA 0xff00

/* header of a block: the bytes which are the same in each block (gzip
   header with the FEXTRA flag, XLEN 6, the BC extra subfield of length 2),
   followed by BSIZE (size of the block minus 1, 16 bits, little endian) */
#define BGZF_HEADER_BYTES \
  0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, \
  'B', 'C', 0x02, 0x00
#define BGZF_BSIZE_OFFSET 16
#define BGZF_HEADER_LEN (BGZF_BSIZE_OFFSET + 2)

/* footer (CRC32, size) of a block */
#define BGZF_FOOTER_LEN 8
#define BGZF_MAX_BLOCK 0x10000

static const uint8_t bgzf_header[BGZF_BSIZE_OFFSET] = {BGZF_HEADER_BYTES};

/* BSIZE 27, empty deflate data, CRC32 0, size 0 */
const uint8_t bgzf_eof[BGZF_EOF_LEN] = {BGZF_HEADER_BYTES,
  0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};

struct bgzf_compressor_t {
  z_stream stream;
  uint8_t block[BGZF_MAX_BLOCK];
};

bgzf_compressor_t* bgzf_compressor_new(int level) {
  bgzf_compressor_t *self;
  XCALLOC(self, 1, sizeof(bgzf_compressor_t));
  /* raw deflate: the gzip header and footer are written here */
  if (deflateInit2(&self->stream, level, Z_DEFLATED, -15, 8,
                   Z_DEFAULT_STRATEGY) != Z_OK)
    fatal_printf("Initialization of the BGZF compression failed\n");
  return self;
}

void bgzf_compressor_delete(bgzf_compressor_t *self) {
  if (self == NULL)
    return;
  deflateEnd(&self->stream);
  free(self);
}

static void put_le32(uint8_t *p, uint32_t v) {
  for (int i = 0; i < 4; i++)
    p[i] = (uint8_t)(v >> (8 * i));
}

/* compress a block of at most BGZF_BLOCK_DATA bytes */
static void compress_block(bgzf_compressor_t *self, const char *data,
    size_t len, strbuf_t *out) {
  uint8_t *b = self->block;
  memcpy(b, bgzf_header, sizeof(bgzf_header));
  deflateReset(&self->stream);
  self->stream.next_in = (Bytef *)data;
  self->stream.avail_in = (uInt)len;
  self->stream.next_out = b + BGZF_HEADER_LEN;
  self->stream.avail_out = BGZF_MAX_BLOCK - BGZF_HEADER_LEN - BGZF_FOOTER_LEN;
  if (deflate(&self->stream, Z_FINISH) != Z_STREAM_END)
    fatal_printf("BGZF compression failed\n");
  size_t block_len = BGZF_HEADER_LEN + self->stream.total_out +
                     BGZF_FOOTER_LEN;
  b[BGZF_BSIZE_OFFSET] = (uint8_t)((block_len - 1) & 0xff);
  b[BGZF_BSIZE_OFFSET + 1] = (uint8_t)((block_len - 1) >> 8);
  uint8_t *footer = b + block_len - BGZF_FOOTER_LEN;
  put_le32(footer, (uint32_t)crc32(crc32(0L, Z_NULL, 0), (const Bytef *)data,
                                   (uInt)len));
  put_le32(footer + 4, (uint32_t)len);
  strbuf_append(out, (const char *)b, block_len);
}

void bgzf_compressor_run(bgzf_compressor_t *self, const char *data,
    size_t len, strbuf_t *out) {
  assert(self != NULL);
  assert(out != NULL);
  for (size_t from = 0; from < len; from += BGZF_BLOCK_DATA)
    compress_block(self, data + from,
        len - from < BGZF_BLOCK_DATA ? len - from : BGZF_BLOCK_DATA, out);
}

/* --- unit tests --- */

void bgzf_compressor_test(const char *tmpdir) {
  /* data larger than a block, with an incompressible part */
  size_t len = 3 * BGZF_BLOCK_DATA + 1000;
  char *data;
  XMALLOC(data, len);
  srand(1);
  for (size_t i = 0; i < len; i++)
    data[i] = i < BGZF_BLOCK_DATA ? (char)(rand() & 0xff) : "ACGT"[i % 4];
  bgzf_compressor_t *c = bgzf_compressor_new(6);
  strbuf_t out;
  strbuf_init(&out);
  bgzf_compressor_run(c, data, len, &out);
  strbuf_append(&out, (const char *)bgzf_eof, BGZF_EOF_LEN);
  /* the blocks: BSIZE gives the size of each block */
  size_t pos = 0, n_blocks = 0;
  while (pos < out.len) {
    const uint8_t *b = (const uint8_t *)out.s + pos;
    test_assert(memcmp(b, bgzf_header, sizeof(bgzf_header)) == 0);
    size_t bsize = (size_t)(b[BGZF_BSIZE_OFFSET] |
                            (b[BGZF_BSIZE_OFFSET + 1] << 8)) + 1;
    test_assert(bsize <= BGZF_MAX_BLOCK);
    pos += bsize;
    n_blocks++;
  }
  test_eq(pos, out.len, size_t, "%zu");
  test_eq(n_blocks, (size_t)5, size_t, "%zu");
  /* the file is a valid gzip file with the original content */
  char *path = strcat_alloc(tmpdir, "/cbaia-tests-bgzf.gz");
  FILE *f = fopen(path, "wb");
  test_assert(f != NULL);
  fwrite(out.s, 1, out.len, f);
  fclose(f);
  gzFile gz = gzopen(path, "rb");
  test_assert(gz != NULL);
  char *back;
  XMALLOC(back, len + 1);
  test_eq(gzread(gz, back, (unsigned)len + 1), (int)len, int, "%d");
  gzclose(gz);
  test_assert(memcmp(back, data, len) == 0);
  unlink(path);
  free(path);
  free(back);
  free(data);
  strbuf_free(&out);
  bgzf_compressor_delete(c);
}
