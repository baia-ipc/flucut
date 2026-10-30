/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "alloc.h"
#include "bam-reader.h"
#include "bgzf-compressor.h"
#include "file.h"
#include "logger.h"
#include "outfiles.h"
#include "reads-reader.h"
#include "strutils.h"
#include "test.h"
#include "test-bam.h"

struct outfiles_t {
  char *prefix;   /* <outdir>/<basename>. */
  FILE **files;
  char **names;   /* name of each destination, once its file is open */
  size_t n_destinations;
  bool append;
  bool bam;            /* false: FASTQ output */
  char **bam_headers;  /* BAM: header of each destination */
  size_t *bam_header_lens;
  Logger *logger;
};

/* create the handler with the given prefix of the file names (copied) */
static outfiles_t* new_with_prefix(const char *prefix,
    size_t n_destinations, Logger *logger) {
  outfiles_t *self;
  XMALLOC(self, sizeof(outfiles_t));
  self->prefix = xstrdup(prefix);
  XCALLOC(self->files, n_destinations + 1, sizeof(FILE *));
  XCALLOC(self->names, n_destinations + 1, sizeof(char *));
  self->n_destinations = n_destinations;
  self->append = false;
  self->bam = false;
  XCALLOC(self->bam_headers, n_destinations + 1, sizeof(char *));
  XCALLOC(self->bam_header_lens, n_destinations + 1, sizeof(size_t));
  self->logger = logger;
  return self;
}

outfiles_t* outfiles_new(const char *outdir, const char *reads_file,
    size_t n_destinations, Logger *logger) {
  assert(outdir != NULL);
  assert(reads_file != NULL);
  xmakedirs(outdir, logger);
  const char *base = strrchr(reads_file, '/');
  base = base == NULL ? reads_file : base + 1;
  size_t len = strlen(base) - reads_reader_filename_suffix_len(base);
  size_t dirlen = strlen(outdir);
  char *prefix;
  XMALLOC(prefix, dirlen + len + 3);
  memcpy(prefix, outdir, dirlen);
  prefix[dirlen] = '/';
  memcpy(prefix + dirlen + 1, base, len);
  prefix[dirlen + 1 + len] = '.';
  prefix[dirlen + 2 + len] = '\0';
  outfiles_t *self = new_with_prefix(prefix, n_destinations, logger);
  free(prefix);
  return self;
}

outfiles_t* outfiles_new_group(const char *outdir, const char *group,
    size_t n_destinations, Logger *logger) {
  assert(outdir != NULL);
  assert(group != NULL);
  char *dir = path_join(outdir, group);
  xmakedirs(dir, logger);
  char *prefix = strcat_alloc(dir, "/");
  outfiles_t *self = new_with_prefix(prefix, n_destinations, logger);
  free(prefix);
  free(dir);
  return self;
}

void outfiles_set_append(outfiles_t *self, bool append) {
  assert(self != NULL);
  self->append = append;
}

void outfiles_set_bam(outfiles_t *self, size_t destination,
    const char *header, size_t len) {
  assert(self != NULL);
  assert(destination < self->n_destinations);
  assert(header != NULL);
  self->bam = true;
  free(self->bam_headers[destination]);
  XMALLOC(self->bam_headers[destination], len + 1);
  memcpy(self->bam_headers[destination], header, len);
  self->bam_header_lens[destination] = len;
}

static void write_or_exit(const outfiles_t *self, FILE *file,
    const char *data, size_t len, const char *name) {
  LCHECK(fwrite(data, 1, len, file) == len, self->logger,
         "Writing the output of %s failed, %s\n", name, strerror(errno));
}

/* is the BAM file non-empty (it has a header); its end-of-file marker is
   removed, so that data can be appended */
static bool prepare_bam_append(const char *path, Logger *logger) {
  struct stat st;
  if (stat(path, &st) != 0 || st.st_size == 0)
    return false;
  if (bam_reader_has_eof_marker(path))
    LCHECK(truncate(path, st.st_size - BGZF_EOF_LEN) == 0, logger,
           "Truncating file '%s' failed, %s\n", path, strerror(errno));
  return true;
}

void outfiles_write(outfiles_t *self, size_t destination, const char *name,
    const char *data, size_t len) {
  assert(self != NULL);
  assert(destination < self->n_destinations);
  if (self->files[destination] == NULL) {
    const char *sfx = self->bam ? OUTFILES_BAM_SFX : OUTFILES_FASTQ_SFX;
    char *base = strcat_alloc(self->prefix, name);
    char *path = strcat_alloc(base, sfx);
    bool has_header = self->bam && self->append &&
                      prepare_bam_append(path, self->logger);
    free(path);
    self->files[destination] = xfopen(base, sfx, self->append ? "ab" : "wb",
                                      "output file", self->logger);
    free(base);
    self->names[destination] = xstrdup(name);
    if (self->bam && !has_header) {
      assert(self->bam_headers[destination] != NULL);
      write_or_exit(self, self->files[destination],
                    self->bam_headers[destination],
                    self->bam_header_lens[destination], name);
    }
  }
  write_or_exit(self, self->files[destination], data, len, name);
}

void outfiles_delete(outfiles_t *self) {
  if (self == NULL)
    return;
  for (size_t i = 0; i < self->n_destinations; i++) {
    if (self->files[i] == NULL)
      continue;
    if (self->bam)
      write_or_exit(self, self->files[i], (const char *)bgzf_eof,
                    BGZF_EOF_LEN, self->names[i]);
    xfclose(&self->files[i], self->logger);
  }
  free(self->files);
  for (size_t i = 0; i < self->n_destinations; i++) {
    free(self->names[i]);
    free(self->bam_headers[i]);
  }
  free(self->names);
  free(self->bam_headers);
  free(self->bam_header_lens);
  free(self->prefix);
  free(self);
}

/* --- unit tests --- */

/* content of a file (allocated), and its length */
static char *slurp(const char *path, size_t *len) {
  FILE *f = fopen(path, "rb");
  test_assert(f != NULL);
  char *buf;
  XMALLOC(buf, 64);
  *len = fread(buf, 1, 64, f);
  fclose(f);
  return buf;
}

#define TEST_BAM_TEXT "@HD\tVN:1.6\tSO:unknown\n"


/* check a BAM file: its header, the names of its reads (separated by
   spaces) and a single end-of-file marker, at the end */
static void check_bam(const char *path, const char *names) {
  FILE *f = fopen(path, "rb");
  test_assert(f != NULL);
  char *content;
  size_t len;
  read_file(f, &content, &len, NULL);
  fclose(f);
  size_t n_eof = 0;
  for (size_t i = 0; i + BGZF_EOF_LEN <= len; i++)
    if (memcmp(content + i, bgzf_eof, BGZF_EOF_LEN) == 0)
      n_eof++;
  test_eq_size(n_eof, 1);
  test_assert(bam_reader_has_eof_marker(path));
  free(content);
  bam_reader_t *r = bam_reader_new(path, NULL);
  test_eq_str(bam_reader_header(r), TEST_BAM_TEXT);
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  strbuf_t found;
  strbuf_init(&found);
  while (bam_reader_read_batch(r, &batch, 2) > 0)
    for (size_t i = 0; i < batch.n_records; i++) {
      if (found.len > 0)
        strbuf_append_char(&found, ' ');
      strbuf_append(&found, batch.data.s + batch.records[i].title_from,
                    batch.records[i].title_len);
    }
  strbuf_append_char(&found, '\0');
  test_eq_str(found.s, names);
  strbuf_free(&found);
  fastq_batch_free(&batch);
  bam_reader_delete(r);
}

/* write the read <title> to the destination <dest> (name <name>) of the
   BAM output files of the input data/run3.bam (2 destinations) */
static void write_bam(const char *outdir, bool append, size_t dest,
    const char *name, const char *title, Logger *logger) {
  strbuf_t data;
  strbuf_init(&data);
  outfiles_t *o = outfiles_new(outdir, "data/run3.bam", 2, logger);
  outfiles_set_append(o, append);
  test_bam_data(&data, TEST_BAM_TEXT, NULL, 0);
  outfiles_set_bam(o, 0, data.s, data.len);
  outfiles_set_bam(o, 1, data.s, data.len);
  strbuf_clear(&data);
  test_bam_data(&data, NULL, &title, 1);
  outfiles_write(o, dest, name, data.s, data.len);
  outfiles_delete(o);
  strbuf_free(&data);
}

static void test_bam(const char *outdir, Logger *logger) {
  /* the suffix .bam of the input is removed; the files are created when
     written, with the header; the end-of-file marker is written when
     they are closed */
  write_bam(outdir, false, 0, "S1", "r1", logger);
  char *path1 = strcat_alloc(outdir, "/run3.S1" OUTFILES_BAM_SFX);
  char *path2 = strcat_alloc(outdir, "/run3.S2" OUTFILES_BAM_SFX);
  check_bam(path1, "r1");
  test_assert(!file_exists(path2, ""));
  /* append: the end-of-file marker is removed, the header is not
     written again; a new file gets the header */
  write_bam(outdir, true, 0, "S1", "r2", logger);
  check_bam(path1, "r1 r2");
  write_bam(outdir, true, 1, "S2", "r3", logger);
  check_bam(path2, "r3");
  /* append to a file without the end-of-file marker (e.g. written by an
     interrupted run): nothing is removed */
  struct stat st;
  test_assert(stat(path2, &st) == 0);
  test_assert(truncate(path2, st.st_size - BGZF_EOF_LEN) == 0);
  write_bam(outdir, true, 1, "S2", "r4", logger);
  check_bam(path2, "r3 r4");
  /* append to an empty file: the header is written */
  FILE *f = fopen(path2, "wb");
  test_assert(f != NULL);
  fclose(f);
  write_bam(outdir, true, 1, "S2", "r5", logger);
  check_bam(path2, "r5");
  /* without append, an existing file is overwritten */
  write_bam(outdir, false, 0, "S1", "r6", logger);
  check_bam(path1, "r6");
  unlink(path1);
  unlink(path2);
  free(path1);
  free(path2);
}

void outfiles_test(const char *tmpdir) {
  Logger *logger = new_logger(stderr, LOG_WARN);
  char *outdir = strcat_alloc(tmpdir, "/cbaia-tests-outfiles");
  outfiles_t *o = outfiles_new_group(outdir, "g1", 2, logger);
  outfiles_write(o, 1, "S2", "ab", 2);
  outfiles_write(o, 1, "S2", "c", 1);
  outfiles_delete(o);
  char *path = strcat_alloc(outdir, "/g1/S2" OUTFILES_FASTQ_SFX);
  char *nopath = strcat_alloc(outdir, "/g1/S1" OUTFILES_FASTQ_SFX);
  size_t len;
  char *content = slurp(path, &len);
  test_assert(len == 3 && memcmp(content, "abc", 3) == 0);
  free(content);
  test_assert(!file_exists(nopath, ""));
  /* append */
  o = outfiles_new_group(outdir, "g1", 2, logger);
  outfiles_set_append(o, true);
  outfiles_write(o, 1, "S2", "d", 1);
  outfiles_delete(o);
  content = slurp(path, &len);
  test_assert(len == 4 && memcmp(content, "abcd", 4) == 0);
  free(content);
  unlink(path);
  /* input file: the suffix of the input is removed */
  o = outfiles_new(outdir, "data/run1.fastq.gz", 1, logger);
  outfiles_write(o, 0, "S1", "x", 1);
  outfiles_delete(o);
  char *path2 = strcat_alloc(outdir, "/run1.S1" OUTFILES_FASTQ_SFX);
  test_assert(file_exists(path2, ""));
  unlink(path2);
  /* also the other suffixes of the files of reads (e.g. .fq.gz) */
  o = outfiles_new(outdir, "run2.fq.gz", 1, logger);
  outfiles_write(o, 0, "S1", "x", 1);
  outfiles_delete(o);
  char *path3 = strcat_alloc(outdir, "/run2.S1" OUTFILES_FASTQ_SFX);
  test_assert(file_exists(path3, ""));
  unlink(path3);
  free(path3);
  o = outfiles_new(outdir, "run4.fastq", 1, logger);
  outfiles_write(o, 0, "S1", "x", 1);
  outfiles_delete(o);
  path3 = strcat_alloc(outdir, "/run4.S1" OUTFILES_FASTQ_SFX);
  test_assert(file_exists(path3, ""));
  unlink(path3);
  free(path3);
  test_bam(outdir, logger);
  char *groupdir = strcat_alloc(outdir, "/g1");
  rmdir(groupdir);
  rmdir(outdir);
  free(groupdir);
  free(path2);
  free(nopath);
  free(path);
  free(outdir);
  delete_logger(logger);
}
