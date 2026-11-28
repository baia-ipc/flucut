/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include "alloc.h"
#include "bam-format.h"
#include "bam-reader.h"
#include "bgzf-compressor.h"
#include "sequtil.h"
#include "strbuf.h"
#include "strutils.h"
#include "test.h"

/* the header text and the records are read in parts of this size */
#define BAM_READ_CHUNK (1 << 20)

/* size of the buffer of zlib for reading the file */
#define BAM_GZBUFFER_SIZE (1 << 17)

/* missing qualities are written as Q0 */
#define QUAL_MISSING_CHAR '!'

struct bam_reader_t {
  gzFile file;
  char *filename;
  uint8_t *rec;    /* current record */
  size_t rec_cap;
  char *rc;        /* reverse complement buffer */
  size_t rc_cap;
  bool done;
  char *header;    /* SAM header text */
  Logger *logger;
};

static uint32_t le32(const uint8_t *p) {
  return (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) |
         ((uint32_t)p[3] << 24);
}

static uint16_t le16(const uint8_t *p) {
  return (uint16_t)(p[0] | (p[1] << 8));
}

static gzFile open_file(const char *filename, Logger *logger) {
  gzFile file = gzopen(filename, "rb");
  LCHECK(file != NULL, logger, "Opening file '%s' failed, %s\n", filename,
         strerror(errno));
  return file;
}

__attribute__ ((noreturn))
static void format_error(const bam_reader_t *self, const char *msg) {
  lprintf(self->logger, LOG_ERROR, "Invalid BAM file %s: %s\n",
          self->filename, msg);
  exit(EXIT_FAILURE);
}

/* read exactly <len> bytes; false at the end of the file (if no byte was
   read and <eof_ok>), otherwise exits on a truncated file */
static bool read_exact(bam_reader_t *self, void *buf, size_t len,
    bool eof_ok) {
  size_t got = 0;
  while (got < len) {
    int n = gzread(self->file, (char *)buf + got, (unsigned)(len - got));
    if (n < 0) {
      int errnum;
      const char *msg = gzerror(self->file, &errnum);
      format_error(self, errnum == Z_ERRNO ? strerror(errno) : msg);
    }
    if (n == 0) {
      if (got == 0 && eof_ok)
        return false;
      format_error(self, "unexpected end of file");
    }
    got += (size_t)n;
  }
  return true;
}

static void skip_bytes(bam_reader_t *self, size_t len) {
  uint8_t buf[4096];
  while (len > 0) {
    size_t n = len < sizeof(buf) ? len : sizeof(buf);
    read_exact(self, buf, n, false);
    len -= n;
  }
}

/* read exactly <len> bytes into <buf> (of capacity *cap, reallocated as
   needed, also with <extra> bytes after the data), in parts of
   BAM_READ_CHUNK bytes, so that the memory is allocated only for the data
   actually present (the length, from the file, may be wrong); returns the
   buffer */
static void *read_chunked(bam_reader_t *self, void *buf, size_t *cap,
    size_t len, size_t extra) {
  size_t got = 0;
  while (got < len) {
    size_t n = len - got < BAM_READ_CHUNK ? len - got : BAM_READ_CHUNK;
    XCHECKSPACEINCPERC(buf, cap, got + n + extra, 50);
    read_exact(self, (char *)buf + got, n, false);
    got += n;
  }
  XCHECKSPACEINCPERC(buf, cap, len + extra, 0);
  return buf;
}

static uint32_t read_u32(bam_reader_t *self) {
  uint8_t b[4];
  read_exact(self, b, 4, false);
  return le32(b);
}

/* read the header (magic, text, references); returns the text
   (0-terminated) */
static char *read_header(bam_reader_t *self) {
  char magic[BAM_MAGIC_LEN];
  if (!read_exact(self, magic, BAM_MAGIC_LEN, true) ||
      memcmp(magic, BAM_MAGIC, BAM_MAGIC_LEN) != 0)
    format_error(self, "wrong magic number");
  size_t l_text = read_u32(self), cap = 0;
  char *text = read_chunked(self, NULL, &cap, l_text, 1);
  /* the text can be padded with NUL characters */
  text[l_text] = '\0';
  uint32_t n_ref = read_u32(self);
  for (uint32_t i = 0; i < n_ref; i++) {
    skip_bytes(self, read_u32(self));
    read_u32(self);
  }
  return text;
}

bool bam_reader_is_bam(const char *filename, Logger *logger) {
  assert(filename != NULL);
  gzFile file = open_file(filename, logger);
  char magic[BAM_MAGIC_LEN];
  int n = gzread(file, magic, BAM_MAGIC_LEN);
  gzclose(file);
  return n == BAM_MAGIC_LEN &&
         memcmp(magic, BAM_MAGIC, BAM_MAGIC_LEN) == 0;
}

bam_reader_t* bam_reader_new(const char *filename, Logger *logger) {
  assert(filename != NULL);
  bam_reader_t *self;
  XMALLOC(self, sizeof(bam_reader_t));
  self->file = open_file(filename, logger);
  self->logger = logger;
  gzbuffer(self->file, BAM_GZBUFFER_SIZE);
  self->filename = xstrdup(filename);
  self->rec = NULL;
  self->rec_cap = 0;
  self->rc = NULL;
  self->rc_cap = 0;
  self->done = false;
  self->header = read_header(self);
  return self;
}

void bam_reader_delete(bam_reader_t *self) {
  if (self == NULL)
    return;
  gzclose(self->file);
  free(self->filename);
  free(self->rec);
  free(self->rc);
  free(self->header);
  free(self);
}

/* read the next record into self->rec; returns its length (0: end) */
static size_t next_record(bam_reader_t *self) {
  uint8_t b[4];
  if (self->done || !read_exact(self, b, 4, true)) {
    self->done = true;
    return 0;
  }
  size_t len = le32(b);
  if (len < BAM_CORE_LEN)
    format_error(self, "record too short");
  /* (block_size, up to 4 GiB, may be wrong) */
  self->rec = read_chunked(self, self->rec, &self->rec_cap, len, 0);
  return len;
}

/* size of the value of a tag of the given type (0: unknown type) */
static size_t tag_value_size(char type) {
  switch (type) {
    case 'A': case 'c': case 'C': return 1;
    case 's': case 'S': return 2;
    case 'i': case 'I': case 'f': return 4;
    default: return 0;
  }
}

/* tags of the base modifications, which are not kept */
static bool skipped_tag(const uint8_t *p) {
  return p[0] == 'M' && (p[1] == 'M' || p[1] == 'L' || p[1] == 'N');
}

/* append the tags (SAM text format) to <out>, each preceded by a tab */
static void append_tags(const bam_reader_t *self, const uint8_t *p,
    const uint8_t *end, strbuf_t *out) {
  while (p < end) {
    if (end - p < 3)
      format_error(self, "truncated tag");
    char type = (char)p[2];
    const uint8_t *v = p + 3;
    size_t size;
    if (type == 'Z' || type == 'H') {
      const uint8_t *z = memchr(v, '\0', (size_t)(end - v));
      if (z == NULL)
        format_error(self, "unterminated string tag");
      size = (size_t)(z - v) + 1;
    } else if (type == 'B') {
      if (end - v < 5)
        format_error(self, "truncated array tag");
      size_t elsize = tag_value_size((char)v[0]);
      if (elsize == 0 || v[0] == 'A')
        format_error(self, "unknown array type");
      size = 5 + elsize * le32(v + 1);
    } else {
      size = tag_value_size(type);
      if (size == 0)
        format_error(self, "unknown tag type");
    }
    if ((size_t)(end - v) < size)
      format_error(self, "truncated tag");
    if (type != 'B' && !skipped_tag(p)) {
      strbuf_append_char(out, '\t');
      strbuf_append(out, (const char *)p, 2);
      strbuf_append_char(out, ':');
      strbuf_append_char(out, type == 'A' || type == 'Z' || type == 'H' ||
                              type == 'f' ? type : 'i');
      strbuf_append_char(out, ':');
      char num[32];
      switch (type) {
        case 'A': strbuf_append_char(out, (char)v[0]); break;
        case 'Z': case 'H':
          strbuf_append(out, (const char *)v, size - 1); break;
        case 'c': strbuf_append_long(out, (int8_t)v[0]); break;
        case 'C': strbuf_append_long(out, v[0]); break;
        case 's': strbuf_append_long(out, (int16_t)le16(v)); break;
        case 'S': strbuf_append_long(out, le16(v)); break;
        case 'i': strbuf_append_long(out, (int32_t)le32(v)); break;
        case 'I':
          snprintf(num, sizeof(num), "%" PRIu32, le32(v));
          strbuf_append_str(out, num);
          break;
        case 'f': {
          uint32_t u = le32(v);
          float f;
          memcpy(&f, &u, sizeof(f));
          snprintf(num, sizeof(num), "%g", (double)f);
          strbuf_append_str(out, num);
          break;
        }
        default:
          /* not reached: the type was checked by tag_value_size() */
          break;
      }
    }
    p = v + size;
  }
}

/* append the record in self->rec (length <len>) to the batch, unless it is
   a secondary or supplementary alignment; returns true if appended */
static bool append_record(bam_reader_t *self, size_t len,
    fastq_batch_t *batch) {
  const uint8_t *r = self->rec, *end = r + len;
  size_t l_read_name = r[8];
  size_t n_cigar_op = le16(r + 12);
  uint16_t flag = le16(r + 14);
  size_t l_seq = le32(r + 16);
  /* the lengths are checked against the record length, as offsets (no
     overflow: l_read_name and n_cigar_op are small), before computing the
     pointers */
  size_t seq_from = BAM_CORE_LEN + l_read_name + 4 * n_cigar_op;
  size_t seq_bytes = l_seq / 2 + l_seq % 2;
  if (l_read_name == 0 || seq_from > len || seq_bytes > len - seq_from ||
      l_seq > len - seq_from - seq_bytes ||
      r[BAM_CORE_LEN + l_read_name - 1] != '\0')
    format_error(self, "inconsistent record lengths");
  const uint8_t *name = r + BAM_CORE_LEN;
  const uint8_t *seq = r + seq_from;
  const uint8_t *qual = seq + seq_bytes;
  const uint8_t *tags = qual + l_seq;
  if (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
    return false;
  bool reverse = (flag & BAM_FREVERSE) != 0;
  XCHECKSPACEINCPERC(batch->records, &batch->records_cap,
                     (batch->n_records + 1) * sizeof(fastq_record_t), 50);
  fastq_record_t *rec = batch->records + batch->n_records++;
  strbuf_t *data = &batch->data;
  rec->title_from = data->len;
  strbuf_append(data, (const char *)name, l_read_name - 1);
  append_tags(self, tags, end, data);
  rec->title_len = data->len - rec->title_from;
  rec->seq_from = data->len;
  rec->seq_len = l_seq;
  XCHECKSPACEINCPERC(self->rc, &self->rc_cap, l_seq + 1, 50);
  for (size_t i = 0; i < l_seq; i++)
    self->rc[i] = BAM_BASES[(seq[i / 2] >> (i % 2 == 0 ? 4 : 0)) & 0xF];
  if (reverse) {
    /* (reverse_complement() also writes a terminator) */
    reverse_complement(self->rc, l_seq, strbuf_reserve(data, l_seq + 1));
    strbuf_commit(data, l_seq);
  } else
    strbuf_append(data, self->rc, l_seq);
  rec->qual_from = data->len;
  /* (with no sequence, the qualities buffer can be NULL) */
  rec->qual_missing = l_seq > 0 && qual[0] == BAM_QUAL_MISSING;
  if (l_seq > 0) {
    char *q = strbuf_reserve(data, l_seq);
    for (size_t i = 0; i < l_seq; i++)
      q[i] = rec->qual_missing ? QUAL_MISSING_CHAR :
             bam_qual_char(qual[reverse ? l_seq - 1 - i : i]);
    strbuf_commit(data, l_seq);
  }
  return true;
}

const char *bam_reader_header(const bam_reader_t *self) {
  assert(self != NULL);
  return self->header;
}

size_t bam_reader_read_batch(bam_reader_t *self, fastq_batch_t *batch,
    size_t max_records) {
  assert(self != NULL);
  assert(batch != NULL);
  strbuf_clear(&batch->data);
  batch->n_records = 0;
  while (batch->n_records < max_records) {
    size_t len = next_record(self);
    if (len == 0)
      break;
    append_record(self, len, batch);
  }
  return batch->n_records;
}

bool bam_reader_has_eof_marker(const char *filename) {
  assert(filename != NULL);
  FILE *f = fopen(filename, "rb");
  if (f == NULL)
    return false;
  uint8_t buf[BGZF_EOF_LEN];
  bool found = fseek(f, -(long)BGZF_EOF_LEN, SEEK_END) == 0 &&
               fread(buf, 1, BGZF_EOF_LEN, f) == BGZF_EOF_LEN &&
               memcmp(buf, bgzf_eof, BGZF_EOF_LEN) == 0;
  fclose(f);
  return found;
}

size_t bam_reader_count(const char *filename, Logger *logger) {
  assert(filename != NULL);
  bam_reader_t *self = bam_reader_new(filename, logger);
  size_t n = 0;
  while (next_record(self) > 0)
    if ((le16(self->rec + 14) & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) == 0)
      n++;
  bam_reader_delete(self);
  return n;
}

/* --- unit tests --- */

/* append a record (sequence of ACGT, qualities as Phred+33, tags as
   binary data) */
static void put_record(strbuf_t *b, const char *name, uint16_t flag,
    const char *seq, const char *qual, const char *tags, size_t tags_len) {
  size_t l_name = strlen(name) + 1, l_seq = strlen(seq);
  size_t n_cigar = (flag & BAM_FUNMAP) ? 0 : 1;
  bam_put32(b, (uint32_t)(BAM_CORE_LEN + l_name + 4 * n_cigar +
                          (l_seq + 1) / 2 + l_seq + tags_len));
  bam_put32(b, (flag & BAM_FUNMAP) ? 0xFFFFFFFF : 0);
  bam_put32(b, (flag & BAM_FUNMAP) ? 0xFFFFFFFF : 100);
  strbuf_append_char(b, (char)l_name);
  strbuf_append_char(b, (char)((flag & BAM_FUNMAP) ? 0 : 60));
  bam_put16(b, BAM_UNMAPPED_BIN);
  bam_put16(b, (uint16_t)n_cigar);
  bam_put16(b, flag);
  bam_put32(b, (uint32_t)l_seq);
  bam_put32(b, 0xFFFFFFFF);
  bam_put32(b, 0xFFFFFFFF);
  bam_put32(b, 0);
  strbuf_append(b, name, l_name);
  if (n_cigar > 0)
    bam_put32(b, (uint32_t)(l_seq << 4));
  for (size_t i = 0; i < l_seq; i += 2) {
    uint8_t hi = (uint8_t)(strchr(BAM_BASES, seq[i]) - BAM_BASES),
            lo = i + 1 < l_seq ?
                 (uint8_t)(strchr(BAM_BASES, seq[i + 1]) - BAM_BASES) : 0;
    strbuf_append_char(b, (char)((hi << 4) | lo));
  }
  for (size_t i = 0; i < l_seq; i++)
    strbuf_append_char(b, qual == NULL ? (char)BAM_QUAL_MISSING :
                       (char)(qual[i] - 33));
  strbuf_append(b, tags, tags_len);
}

static void check_record(const fastq_batch_t *batch, size_t num,
    const char *title, const char *seq, const char *qual) {
  const fastq_record_t *r = batch->records + num;
  const char *d = batch->data.s;
  test_eq(r->title_len, strlen(title), size_t, "%zu");
  test_assert(memcmp(d + r->title_from, title, strlen(title)) == 0);
  test_eq(r->seq_len, strlen(seq), size_t, "%zu");
  test_assert(memcmp(d + r->seq_from, seq, strlen(seq)) == 0);
  test_assert(memcmp(d + r->qual_from, qual, strlen(qual)) == 0);
}

void bam_reader_test(const char *tmpdir) {
  strbuf_t b;
  strbuf_init(&b);
  strbuf_append(&b, BAM_MAGIC, BAM_MAGIC_LEN);
  const char *text = "@HD\tVN:1.6\tSO:unknown\n";
  bam_put32(&b, (uint32_t)strlen(text));
  strbuf_append_str(&b, text);
  bam_put32(&b, 1);
  bam_put32(&b, 5);
  strbuf_append(&b, "chr1", 5);
  bam_put32(&b, 1000);
  /* tags: qs:i (as C), st:Z, du:f (1.5), ch:s (-3), MM:Z, ML:B:C, RG:Z */
  const char tags[] = "qsC\x0c" "stZ2026-01-01\0" "duf\x00\x00\xc0\x3f"
                      "chs\xfd\xff" "MMZC+m?,0;\0" "MLBC\x01\x00\x00\x00\xc8"
                      "RGZrun1\0";
  put_record(&b, "r1", BAM_FUNMAP, "ACGTN", "ABCDE", tags,
             sizeof(tags) - 1);
  put_record(&b, "r2", BAM_FREVERSE, "AACGT", "ABCDE", NULL, 0);
  put_record(&b, "r3", BAM_FSECONDARY, "AAAA", "AAAA", NULL, 0);
  put_record(&b, "r4", BAM_FSUPPLEMENTARY, "AAAA", "AAAA", NULL, 0);
  put_record(&b, "r5", BAM_FUNMAP, "ACG", NULL, NULL, 0);
  /* quality values above 93 (94, 200), clamped to '~' */
  put_record(&b, "r6", BAM_FUNMAP, "ACG", "\x7f\xe9I", NULL, 0);
  /* a record longer than the parts in which it is read */
  size_t long_len = BAM_READ_CHUNK + 5;
  char *long_seq, *long_qual;
  XMALLOC(long_seq, long_len + 1);
  XMALLOC(long_qual, long_len + 1);
  memset(long_seq, 'C', long_len);
  memset(long_qual, 'I', long_len);
  long_seq[long_len] = long_qual[long_len] = '\0';
  put_record(&b, "r7", BAM_FUNMAP, long_seq, long_qual, NULL, 0);
  char *path = strcat_alloc(tmpdir, "/cbaia-tests-reads.bam");
  gzFile f = gzopen(path, "wb");
  test_assert(f != NULL);
  test_eq(gzwrite(f, b.s, (unsigned)b.len), (int)b.len, int, "%d");
  gzclose(f);
  strbuf_free(&b);
  test_assert(bam_reader_is_bam(path, NULL));
  bam_reader_t *h = bam_reader_new(path, NULL);
  test_assert(strcmp(bam_reader_header(h), text) == 0);
  bam_reader_delete(h);
  test_assert(!bam_reader_has_eof_marker(path));
  FILE *raw = fopen(path, "ab");
  test_assert(raw != NULL);
  fwrite(bgzf_eof, 1, BGZF_EOF_LEN, raw);
  fclose(raw);
  test_assert(bam_reader_has_eof_marker(path));
  test_eq(bam_reader_count(path, NULL), (size_t)5, size_t, "%zu");
  bam_reader_t *r = bam_reader_new(path, NULL);
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  test_eq(bam_reader_read_batch(r, &batch, 2), (size_t)2, size_t, "%zu");
  check_record(&batch, 0, "r1\tqs:i:12\tst:Z:2026-01-01\tdu:f:1.5\t"
               "ch:i:-3\tRG:Z:run1", "ACGTN", "ABCDE");
  check_record(&batch, 1, "r2", "ACGTT", "EDCBA");
  test_assert(!batch.records[0].qual_missing);
  test_eq(bam_reader_read_batch(r, &batch, 2), (size_t)2, size_t, "%zu");
  check_record(&batch, 0, "r5", "ACG", "!!!");
  test_assert(batch.records[0].qual_missing);
  check_record(&batch, 1, "r6", "ACG", "~~I");
  test_assert(!batch.records[1].qual_missing);
  test_eq(bam_reader_read_batch(r, &batch, 2), (size_t)1, size_t, "%zu");
  check_record(&batch, 0, "r7", long_seq, long_qual);
  test_eq(bam_reader_read_batch(r, &batch, 2), (size_t)0, size_t, "%zu");
  free(long_seq);
  free(long_qual);
  fastq_batch_free(&batch);
  bam_reader_delete(r);
  /* a FASTQ file is not a BAM file */
  f = gzopen(path, "wb");
  gzputs(f, "@r1\nACGT\n+\nAAAA\n");
  gzclose(f);
  test_assert(!bam_reader_is_bam(path, NULL));
  unlink(path);
  free(path);
}
