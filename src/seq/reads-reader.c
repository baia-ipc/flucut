/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include "alloc.h"
#include "bam-reader.h"
#include "file.h"
#include "reads-reader.h"
#include "strutils.h"
#include "test.h"
#include "test-bam.h"

#define READS_READER_BUFSIZE (1 << 20)

/* first two bytes of a gzip file */
#define GZIP_MAGIC_1 0x1f
#define GZIP_MAGIC_2 0x8b

struct reads_reader_t {
  fastq_reader_t *fastq;
  bam_reader_t *bam;
};

reads_reader_t* reads_reader_new(const char *filename, Logger *logger) {
  assert(filename != NULL);
  reads_reader_t *self;
  XMALLOC(self, sizeof(reads_reader_t));
  self->fastq = NULL;
  self->bam = NULL;
  if (bam_reader_is_bam(filename, logger))
    self->bam = bam_reader_new(filename, logger);
  else
    self->fastq = fastq_reader_new(filename, logger);
  return self;
}

bool reads_reader_is_bam(const reads_reader_t *self) {
  assert(self != NULL);
  return self->bam != NULL;
}

size_t reads_reader_read_batch(reads_reader_t *self, fastq_batch_t *batch,
    size_t max_records) {
  assert(self != NULL);
  return self->bam != NULL ?
         bam_reader_read_batch(self->bam, batch, max_records) :
         fastq_reader_read_batch(self->fastq, batch, max_records);
}

void reads_reader_delete(reads_reader_t *self) {
  if (self == NULL)
    return;
  bam_reader_delete(self->bam);
  fastq_reader_delete(self->fastq);
  free(self);
}

/* the suffixes of the files of reads (also removed from the input name by
   outfiles_new(), thus a single list) */
static const char *const reads_filename_suffixes[] =
  {".fastq", ".fastq.gz", ".fq", ".fq.gz", ".bam", NULL};

size_t reads_reader_filename_suffix_len(const char *name) {
  for (const char *const *sfx = reads_filename_suffixes; *sfx != NULL; sfx++)
    if (str_has_suffix(name, *sfx))
      return strlen(*sfx);
  return 0;
}

bool reads_reader_is_reads_filename(const char *name) {
  return reads_reader_filename_suffix_len(name) > 0;
}

size_t reads_reader_count(const char *filename, Logger *logger) {
  return bam_reader_is_bam(filename, logger) ?
         bam_reader_count(filename, logger) :
         fastq_reader_count(filename, logger);
}

/* whether an uncompressed file ends with a newline (or is empty); false
   if it cannot be read; *is_gzip: whether it is gzip-compressed instead
   (and then it is not read further) */
static bool plain_file_complete(const char *path, bool *is_gzip) {
  *is_gzip = false;
  FILE *f = fopen(path, "rb");
  if (f == NULL)
    return false;
  unsigned char magic[2];
  size_t n = fread(magic, 1, 2, f);
  bool complete;
  if (n == 2 && magic[0] == GZIP_MAGIC_1 && magic[1] == GZIP_MAGIC_2) {
    *is_gzip = true;
    complete = false;
  } else if (n == 0)
    complete = !ferror(f);
  else
    complete = fseek(f, -1, SEEK_END) == 0 && fgetc(f) == '\n';
  fclose(f);
  return complete;
}

bool reads_reader_file_complete(const char *path) {
  if (str_has_suffix(path, ".bam"))
    return bam_reader_has_eof_marker(path);
  bool is_gzip;
  bool complete = plain_file_complete(path, &is_gzip);
  if (!is_gzip)
    return complete;
  /* the last character of a gzip file is known only by decompressing the
     whole file (no index; the trailer has only the CRC and the length):
     the time is that of decompressing the file once (as reading it); it
     also verifies the integrity of the data (CRC); thus, a program
     polling the files being written should check their age first, and
     call it once for each file */
  gzFile f = gzopen(path, "rb");
  if (f == NULL)
    return false;
  char *buf;
  XMALLOC(buf, READS_READER_BUFSIZE);
  int n;
  char last = '\n';
  while ((n = gzread(f, buf, READS_READER_BUFSIZE)) > 0)
    last = buf[n - 1];
  int err;
  gzerror(f, &err);
  free(buf);
  gzclose(f);
  return n == 0 && err == Z_OK && last == '\n';
}

/* --- unit tests --- */

static char *write_gz(const char *tmpdir, const char *name, const char *data,
    size_t len) {
  char *path = path_join(tmpdir, name);
  gzFile f = gzopen(path, "wb");
  test_assert(f != NULL);
  test_assert(gzwrite(f, data, (unsigned)len) == (int)len);
  gzclose(f);
  return path;
}

/* write a BAM file with the reads r1 and r2 (with the end-of-file marker
   if <eof>) */
static char *write_bam(const char *tmpdir, const char *name, bool eof) {
  const char *const titles[] = {"r1", "r2"};
  return test_bam_file(tmpdir, name, "@HD\tVN:1.6\tSO:unknown\n", titles,
                       2, eof);
}

static void test_bam(const char *tmpdir, Logger *logger) {
  char *path = write_bam(tmpdir, "cbaia-tests-reads-reader.bam", true);
  test_assert(reads_reader_file_complete(path));
  test_eq_size(reads_reader_count(path, logger), 2);
  reads_reader_t *r = reads_reader_new(path, logger);
  test_assert(reads_reader_is_bam(r));
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  test_eq_size(reads_reader_read_batch(r, &batch, 1), 1);
  test_eq_size(batch.records[0].title_len, 2);
  test_assert(memcmp(batch.data.s + batch.records[0].title_from, "r1", 2)
              == 0);
  test_eq_size(reads_reader_read_batch(r, &batch, 10), 1);
  test_assert(memcmp(batch.data.s + batch.records[0].title_from, "r2", 2)
              == 0);
  test_assert(memcmp(batch.data.s + batch.records[0].seq_from, "ACGT", 4)
              == 0);
  test_eq_size(reads_reader_read_batch(r, &batch, 10), 0);
  fastq_batch_free(&batch);
  reads_reader_delete(r);
  unlink(path);
  free(path);
  /* without the end-of-file marker: not complete (but readable) */
  path = write_bam(tmpdir, "cbaia-tests-reads-reader2.bam", false);
  test_assert(!reads_reader_file_complete(path));
  test_eq_size(reads_reader_count(path, logger), 2);
  unlink(path);
  free(path);
  test_assert(!reads_reader_file_complete("/nonexistent/reads.bam"));
}

void reads_reader_test(const char *tmpdir) {
  Logger *logger = new_logger(stderr, LOG_WARN);
  test_assert(reads_reader_is_reads_filename("a.fastq.gz"));
  test_assert(reads_reader_is_reads_filename("a.fq"));
  test_assert(reads_reader_is_reads_filename("a.bam"));
  test_assert(!reads_reader_is_reads_filename("a.txt"));
  test_assert(!reads_reader_is_reads_filename(".fastq"));
  test_eq_size(reads_reader_filename_suffix_len("a.fq.gz"), 6);
  test_eq_size(reads_reader_filename_suffix_len("a.fastq"), 6);
  test_eq_size(reads_reader_filename_suffix_len("a.gz"), 0);
  const char *fq = "@r1\nACGT\n+\nIIII\n@r2\nGG\n+\nII\n";
  char *path = write_gz(tmpdir, "cbaia-tests-reads-reader.fastq.gz", fq,
                        strlen(fq));
  test_assert(reads_reader_file_complete(path));
  test_eq_size(reads_reader_count(path, logger), 2);
  reads_reader_t *r = reads_reader_new(path, logger);
  test_assert(!reads_reader_is_bam(r));
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  test_eq_size(reads_reader_read_batch(r, &batch, 10), 2);
  test_eq_size(reads_reader_read_batch(r, &batch, 10), 0);
  fastq_batch_free(&batch);
  reads_reader_delete(r);
  unlink(path);
  free(path);
  /* last line incomplete */
  path = write_gz(tmpdir, "cbaia-tests-reads-reader2.fastq.gz", fq,
                  strlen(fq) - 1);
  test_assert(!reads_reader_file_complete(path));
  test_eq_size(reads_reader_count(path, logger), 2);
  unlink(path);
  free(path);
  /* multi-line records, blank lines: the records are counted, not lines */
  const char *fq2 = "@r1\nAC\nGT\n+\nII\nII\n\n@r2\nG\n+\nI\n";
  path = write_gz(tmpdir, "cbaia-tests-reads-reader3.fastq.gz", fq2,
                  strlen(fq2));
  test_eq_size(reads_reader_count(path, logger), 2);
  unlink(path);
  free(path);
  test_assert(!reads_reader_file_complete("/nonexistent/reads.fastq.gz"));
  /* uncompressed: complete if empty or ending with a newline */
  path = path_join(tmpdir, "cbaia-tests-reads-reader.fastq");
  FILE *f = fopen(path, "w");
  test_assert(f != NULL);
  fclose(f);
  test_assert(reads_reader_file_complete(path));
  f = fopen(path, "w");
  test_assert(f != NULL);
  fputs(fq, f);
  fclose(f);
  test_assert(reads_reader_file_complete(path));
  test_assert(truncate(path, (off_t)strlen(fq) - 1) == 0);
  test_assert(!reads_reader_file_complete(path));
  test_assert(truncate(path, 1) == 0);
  test_assert(!reads_reader_file_complete(path));
  unlink(path);
  free(path);
  test_bam(tmpdir, logger);
  delete_logger(logger);
}
