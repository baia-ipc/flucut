/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include "alloc.h"
#include "fastq-reader.h"
#include "strutils.h"
#include "test.h"

#define FASTQ_READER_BUFSIZE (1 << 20)

/* maximal length of the description of a format error */
#define FASTQ_READER_MSG_LEN 256

struct fastq_reader_t {
  char *filename;
  gzFile file;
  char *buf;
  size_t buf_pos, buf_len;
  bool eof;
  size_t lineno;
  /* current line: the title line of the next record, once read */
  strbuf_t line;
  bool started, done;
  Logger *logger;
};

static fastq_reader_t* open_reader(const char *filename, Logger *logger) {
  assert(filename != NULL);
  fastq_reader_t *self;
  XMALLOC(self, sizeof(fastq_reader_t));
  self->filename = xstrdup(filename);
  /* gzopen reads also uncompressed files */
  self->file = gzopen(filename, "rb");
  LCHECK(self->file != NULL, logger, "Opening file '%s' failed, %s\n",
         filename, strerror(errno));
  self->logger = logger;
  gzbuffer(self->file, FASTQ_READER_BUFSIZE);
  XMALLOC(self->buf, FASTQ_READER_BUFSIZE);
  self->buf_pos = 0;
  self->buf_len = 0;
  self->eof = false;
  self->lineno = 0;
  strbuf_init(&self->line);
  self->started = false;
  self->done = false;
  return self;
}

fastq_reader_t* fastq_reader_new(const char *filename, Logger *logger) {
  fastq_reader_t *self = open_reader(filename, logger);
  lprintf(logger, LOG_INFO, "Reading %s\n", filename);
  return self;
}

void fastq_reader_delete(fastq_reader_t *self) {
  if (self == NULL)
    return;
  gzclose(self->file);
  free(self->buf);
  strbuf_free(&self->line);
  free(self->filename);
  free(self);
}

void fastq_batch_init(fastq_batch_t *batch) {
  strbuf_init(&batch->data);
  batch->records = NULL;
  batch->n_records = 0;
  batch->records_cap = 0;
}

void fastq_batch_free(fastq_batch_t *batch) {
  strbuf_free(&batch->data);
  free(batch->records);
  fastq_batch_init(batch);
}

static void fill_buffer(fastq_reader_t *self) {
  int n = gzread(self->file, self->buf, FASTQ_READER_BUFSIZE);
  if (n < 0) {
    int errnum;
    const char *msg = gzerror(self->file, &errnum);
    lprintf(self->logger, LOG_ERROR, "Reading file '%s' failed, %s\n",
            self->filename, errnum == Z_ERRNO ? strerror(errno) : msg);
    exit(EXIT_FAILURE);
  }
  self->buf_pos = 0;
  self->buf_len = (size_t)n;
  if (n == 0)
    self->eof = true;
}

/* append the next line (without newline) to out; returns false at the
   end of the file (if no character was read) */
static bool read_line(fastq_reader_t *self, strbuf_t *out) {
  bool any = false;
  while (true) {
    if (self->buf_pos == self->buf_len) {
      if (self->eof)
        break;
      fill_buffer(self);
      if (self->eof)
        break;
    }
    char *start = self->buf + self->buf_pos;
    size_t avail = self->buf_len - self->buf_pos;
    char *nl = memchr(start, '\n', avail);
    any = true;
    if (nl != NULL) {
      strbuf_append(out, start, (size_t)(nl - start));
      self->buf_pos += (size_t)(nl - start) + 1;
      self->lineno++;
      return true;
    }
    strbuf_append(out, start, avail);
    self->buf_pos = self->buf_len;
  }
  if (any)
    self->lineno++;
  return any;
}

/* remove trailing whitespace from data[from..] */
static void rstrip(strbuf_t *data, size_t from) {
  size_t len = data->len;
  while (len > from && is_space(data->s[len - 1]))
    len--;
  strbuf_truncate(data, len);
}

__attribute__ ((format (printf, 2, 3), noreturn))
static void format_error(const fastq_reader_t *self, const char *fmt, ...) {
  char msg[FASTQ_READER_MSG_LEN];
  va_list ap;
  va_start(ap, fmt);
  vsnprintf(msg, sizeof(msg), fmt, ap);
  va_end(ap);
  lprintf(self->logger, LOG_ERROR, "Invalid FASTQ file %s, line %zu: %s\n",
          self->filename, self->lineno, msg);
  exit(EXIT_FAILURE);
}

/* read the next line into self->line; false at the end of the file */
static bool next_line(fastq_reader_t *self) {
  strbuf_clear(&self->line);
  return read_line(self, &self->line);
}

static inline char first_char(const strbuf_t *line) {
  return line->len > 0 ? line->s[0] : '\n';
}

/* append the current line from position <from> to the data, removing
   trailing whitespace */
static void append_line(fastq_reader_t *self, size_t from, strbuf_t *data) {
  size_t piece_from = data->len;
  strbuf_append(data, self->line.s + from, self->line.len - from);
  rstrip(data, piece_from);
}

/* parsing as the Biopython FastqGeneralIterator, see fastq-reader.h */
size_t fastq_reader_read_batch(fastq_reader_t *self, fastq_batch_t *batch,
    size_t max_records) {
  assert(self != NULL);
  assert(batch != NULL);
  strbuf_clear(&batch->data);
  batch->n_records = 0;
  strbuf_t *data = &batch->data;
  if (!self->started) {
    self->started = true;
    if (!next_line(self))
      self->done = true;
  }
  while (!self->done && batch->n_records < max_records) {
    if (first_char(&self->line) != '@')
      format_error(self, "records in FASTQ files should start with '@'");
    XCHECKSPACEINCPERC(batch->records, &batch->records_cap,
                       (batch->n_records + 1) * sizeof(fastq_record_t), 50);
    fastq_record_t *rec = batch->records + batch->n_records;
    /* title */
    rec->title_from = data->len;
    append_line(self, 1, data);
    rec->title_len = data->len - rec->title_from;
    /* sequence, until the "+" line */
    rec->seq_from = data->len;
    while (true) {
      if (!next_line(self))
        format_error(self, data->len > rec->seq_from ?
                     "end of file without quality information" :
                     "unexpected end of file");
      if (first_char(&self->line) == '+')
        break;
      size_t piece_from = data->len;
      append_line(self, 0, data);
      for (size_t i = piece_from; i < data->len; i++)
        if (data->s[i] == ' ' || data->s[i] == '\t')
          format_error(self, "whitespace is not allowed in the sequence");
    }
    rec->seq_len = data->len - rec->seq_from;
    /* the title of the "+" line is optional, but must match */
    size_t plus_title = data->len;
    append_line(self, 1, data);
    if (data->len > plus_title &&
        (data->len - plus_title != rec->title_len ||
         memcmp(data->s + plus_title, data->s + rec->title_from,
                rec->title_len) != 0))
      format_error(self, "sequence and quality captions differ");
    strbuf_truncate(data, plus_title);
    /* qualities, until the title of the next record or the end of file */
    rec->qual_from = data->len;
    rec->qual_missing = false;
    bool any = false;
    while (true) {
      if (!next_line(self)) {
        if (!any)
          format_error(self, "unexpected end of file");
        self->done = true;
        break;
      }
      any = true;
      if (first_char(&self->line) == '@' &&
          data->len - rec->qual_from >= rec->seq_len)
        break;
      append_line(self, 0, data);
    }
    if (data->len - rec->qual_from != rec->seq_len)
      format_error(self, "lengths of sequence and quality values differ "
                   "(%zu and %zu)", rec->seq_len,
                   data->len - rec->qual_from);
    batch->n_records++;
  }
  return batch->n_records;
}

/* the records are parsed (not only the lines counted), so that the count is
   that of the reads read, also for multi-line records, blank lines or a
   missing final newline; the time is anyway mostly that of decompressing */
size_t fastq_reader_count(const char *filename, Logger *logger) {
  fastq_reader_t *self = open_reader(filename, logger);
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  size_t n = 0, n_batch;
  while ((n_batch = fastq_reader_read_batch(self, &batch, 1024)) > 0)
    n += n_batch;
  fastq_batch_free(&batch);
  fastq_reader_delete(self);
  return n;
}

/* --- unit tests --- */

void fastq_reader_test(const char *tmpdir) {
  char *path = strcat_alloc(tmpdir, "/cbaia-tests-reads.fastq.gz");
  gzFile f = gzopen(path, "wb");
  test_assert(f != NULL);
  gzputs(f, "@r1 a\tb  \nACGT\n+\nIIII\n\n@r2\nAC\n+r2\n#$\n@r3\nA\n+\nI");
  gzclose(f);
  fastq_reader_t *reader = fastq_reader_new(path, NULL);
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  test_eq_size(fastq_reader_read_batch(reader, &batch, 2), 2);
  const fastq_record_t *r = batch.records;
  test_assert(strncmp(batch.data.s + r->title_from, "r1 a\tb", 6) == 0);
  test_eq_size(r->title_len, 6);
  test_assert(strncmp(batch.data.s + r->seq_from, "ACGT", 4) == 0);
  test_eq_size(r->seq_len, 4);
  test_assert(strncmp(batch.data.s + r->qual_from, "IIII", 4) == 0);
  test_assert(!r->qual_missing);
  r = batch.records + 1;
  test_eq_size(r->title_len, 2);
  test_assert(strncmp(batch.data.s + r->qual_from, "#$", 2) == 0);
  test_eq_size(fastq_reader_read_batch(reader, &batch, 2), 1);
  test_eq_size(batch.records[0].seq_len, 1);
  test_eq_size(fastq_reader_read_batch(reader, &batch, 2), 0);
  fastq_reader_delete(reader);
  /* multi-line records; quality lines starting with '@' */
  f = gzopen(path, "wb");
  test_assert(f != NULL);
  gzputs(f, "@r1\nACG\nTA\n+\n@I\nIII\n@r2\nAC\n+\n@@\n@r3\nA\n+\nI\n");
  gzclose(f);
  reader = fastq_reader_new(path, NULL);
  test_eq_size(fastq_reader_read_batch(reader, &batch, 10), 3);
  r = batch.records;
  test_eq_size(r->seq_len, 5);
  test_assert(strncmp(batch.data.s + r->seq_from, "ACGTA", 5) == 0);
  test_assert(strncmp(batch.data.s + r->qual_from, "@IIII", 5) == 0);
  r = batch.records + 1;
  test_assert(strncmp(batch.data.s + r->title_from, "r2", 2) == 0);
  test_assert(strncmp(batch.data.s + r->qual_from, "@@", 2) == 0);
  test_eq_size(fastq_reader_read_batch(reader, &batch, 10), 0);
  fastq_reader_delete(reader);
  test_eq_size(fastq_reader_count(path, NULL), 3);
  unlink(path);
  free(path);
  /* uncompressed, with CRLF line ends; records with an empty sequence */
  path = strcat_alloc(tmpdir, "/cbaia-tests-reads.fastq");
  FILE *pf = fopen(path, "w");
  test_assert(pf != NULL);
  fputs("@r1 x\r\nACGT\r\n+r1 x\r\nIIII\r\n@r2\r\n+\r\n"
        "@r3\r\n\r\n+\r\n\r\n@r4\r\nAC\r\n+\r\nII\r\n", pf);
  fclose(pf);
  reader = fastq_reader_new(path, NULL);
  test_eq_size(fastq_reader_read_batch(reader, &batch, 10), 4);
  r = batch.records;
  test_eq_mem(batch.data.s + r->title_from, r->title_len, "r1 x");
  test_eq_mem(batch.data.s + r->seq_from, r->seq_len, "ACGT");
  test_eq_mem(batch.data.s + r->qual_from, r->seq_len, "IIII");
  for (size_t i = 1; i <= 2; i++) {
    r = batch.records + i;
    test_eq_size(r->title_len, 2);
    test_eq_size(r->seq_len, 0);
  }
  r = batch.records + 3;
  test_eq_mem(batch.data.s + r->title_from, r->title_len, "r4");
  test_eq_mem(batch.data.s + r->seq_from, r->seq_len, "AC");
  test_eq_mem(batch.data.s + r->qual_from, r->seq_len, "II");
  test_assert(memchr(batch.data.s, '\r', batch.data.len) == NULL);
  test_eq_size(fastq_reader_read_batch(reader, &batch, 10), 0);
  fastq_reader_delete(reader);
  fastq_batch_free(&batch);
  unlink(path);
  free(path);
}
