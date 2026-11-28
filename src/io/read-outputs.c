/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "alloc.h"
#include "file.h"
#include "input-groups.h"
#include "read-outputs.h"
#include "strbuf.h"
#include "strutils.h"
#include "test.h"

/* the path <outdir>/<group>.<middle>.fastq.gz */
static char* category_path(const char *outdir, const char *group,
    const char *middle) {
  strbuf_t name;
  strbuf_init(&name);
  strbuf_append_str(&name, group);
  strbuf_append_char(&name, '.');
  strbuf_append_str(&name, middle);
  strbuf_append_char(&name, '\0');
  char *path = input_group_output_path(outdir, name.s);
  strbuf_free(&name);
  return path;
}

/* the names of the categories as {full,internal,...} */
static char* category_names(void) {
  strbuf_t names;
  strbuf_init(&names);
  for (size_t k = 0; k < READ_N_CATEGORIES; k++) {
    strbuf_append_char(&names, k == 0 ? '{' : ',');
    strbuf_append_str(&names, read_category_name(k));
  }
  strbuf_append_char(&names, '}');
  strbuf_append_char(&names, '\0');
  char *result = xstrdup(names.s);
  strbuf_free(&names);
  return result;
}

read_outputs_t* read_outputs_new(const char *outdir, const char *group,
    bool by_category) {
  assert(outdir != NULL);
  assert(group != NULL);
  read_outputs_t *self;
  XCALLOC(self, 1, sizeof(read_outputs_t));
  if (!by_category) {
    self->n = 1;
    self->paths[0] = input_group_output_path(outdir, group);
    self->label = xstrdup(self->paths[0]);
    return self;
  }
  self->n = READ_N_CATEGORIES;
  for (size_t k = 0; k < self->n; k++)
    self->paths[k] = category_path(outdir, group, read_category_name(k));
  char *names = category_names();
  self->label = category_path(outdir, group, names);
  free(names);
  return self;
}

void read_outputs_open(read_outputs_t *self, bool append, Logger *logger) {
  assert(self != NULL);
  for (size_t k = 0; k < self->n; k++) {
    assert(self->files[k] == NULL);
    self->files[k] = xfopen(self->paths[k], NULL, append ? "ab" : "wb",
                            "trimmed reads", logger);
  }
}

void read_outputs_write(read_outputs_t *self, size_t k, const char *data,
    size_t len, Logger *logger) {
  assert(self != NULL);
  assert(k < self->n);
  assert(self->files[k] != NULL);
  if (len == 0)
    return;
  LCHECK(fwrite(data, 1, len, self->files[k]) == len, logger,
         "%s: write error\n", self->paths[k]);
}

/* a gzip member without data: header, an empty final block, the CRC and
   the length (both 0) */
static const unsigned char empty_gzip[] = {
    0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0};

void read_outputs_close(read_outputs_t *self, Logger *logger) {
  assert(self != NULL);
  for (size_t k = 0; k < self->n; k++) {
    FILE *f = self->files[k];
    if (f == NULL)
      continue;
    /* (a file without any byte is not a gzip file) */
    LCHECK(fseek(f, 0, SEEK_END) == 0 && (ftell(f) != 0 ||
           fwrite(empty_gzip, 1, sizeof(empty_gzip), f) ==
           sizeof(empty_gzip)), logger, "%s: write error\n",
           self->paths[k]);
    xfclose(&self->files[k], logger);
  }
}

void read_outputs_remove(const read_outputs_t *self, Logger *logger) {
  assert(self != NULL);
  for (size_t k = 0; k < self->n; k++)
    LCHECK(unlink(self->paths[k]) == 0 || errno == ENOENT, logger,
           "Deleting file '%s' failed, %s\n", self->paths[k],
           strerror(errno));
}

bool read_outputs_has_path(const read_outputs_t *self, const char *path) {
  assert(self != NULL);
  for (size_t k = 0; k < self->n; k++)
    if (strcmp(self->paths[k], path) == 0)
      return true;
  return false;
}

void read_outputs_delete(read_outputs_t *self) {
  if (self == NULL)
    return;
  for (size_t k = 0; k < self->n; k++) {
    assert(self->files[k] == NULL);
    free(self->paths[k]);
  }
  free(self->label);
  free(self);
}

/* --- unit tests --- */

/* the size of the file <path>, -1 if it does not exist */
static long test_size(const char *path) {
  FILE *f = fopen(path, "rb");
  if (f == NULL)
    return -1;
  fseek(f, 0, SEEK_END);
  long size = ftell(f);
  fclose(f);
  return size;
}

void read_outputs_test(const char *tmpdir) {
  Logger *logger = new_logger(stderr, LOG_ERROR);
  /* one file */
  read_outputs_t *o = read_outputs_new("out", "bc01", false);
  test_eq_size(o->n, 1);
  test_eq_str(o->paths[0], "out/bc01.fastq.gz");
  test_eq_str(o->label, "out/bc01.fastq.gz");
  test_assert(read_outputs_has_path(o, "out/bc01.fastq.gz"));
  test_assert(!read_outputs_has_path(o, "out/bc01.full.fastq.gz"));
  read_outputs_delete(o);
  /* a file for each category: all created, written, appended to and
     removed */
  o = read_outputs_new(tmpdir, "bc01", true);
  test_eq_size(o->n, READ_N_CATEGORIES);
  test_assert(strstr(o->paths[READ_CATEGORY_FULL],
                     "/bc01.full.fastq.gz") != NULL);
  test_assert(strstr(o->paths[READ_CATEGORY_UNPLACED],
                     "/bc01.unplaced.fastq.gz") != NULL);
  test_assert(strstr(o->paths[READ_CATEGORY_UNKNOWN],
                     "/bc01.unknown.fastq.gz") != NULL);
  test_assert(strstr(o->label,
      "/bc01.{full,internal,deletion,fragment,unplaced,unknown}.fastq.gz")
      != NULL);
  test_assert(read_outputs_has_path(o, o->paths[READ_CATEGORY_DELETION]));
  read_outputs_open(o, false, logger);
  read_outputs_write(o, READ_CATEGORY_INTERNAL, "abc", 3, logger);
  read_outputs_close(o, logger);
  test_eq_long(test_size(o->paths[READ_CATEGORY_INTERNAL]), 3);
  /* (a file without reads: an empty gzip member, written once) */
  test_eq_long(test_size(o->paths[READ_CATEGORY_FULL]), 20);
  read_outputs_open(o, true, logger);
  read_outputs_write(o, READ_CATEGORY_INTERNAL, "de", 2, logger);
  read_outputs_close(o, logger);
  test_eq_long(test_size(o->paths[READ_CATEGORY_INTERNAL]), 5);
  test_eq_long(test_size(o->paths[READ_CATEGORY_FULL]), 20);
  read_outputs_remove(o, logger);
  test_eq_long(test_size(o->paths[READ_CATEGORY_INTERNAL]), -1);
  read_outputs_remove(o, logger);   /* (already removed) */
  read_outputs_delete(o);
  delete_logger(logger);
}
