/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "fasta.h"
#include "file.h"
#include "logger.h"
#include "primers.h"
#include "strbuf.h"
#include "strutils.h"
#include "test.h"

/* append a record (name, description and sequence) to <self> */
static void add_record(fasta_t *self, size_t *cap, const char *name,
    const char *description, const strbuf_t *seq, const char *path,
    Logger *logger) {
  LCHECK(seq->len > 0, logger,
         "%s: the sequence of \"%s\" is empty\n", path, name);
  if (self->n_records == *cap) {
    *cap = *cap == 0 ? 8 : *cap * 2;
    XREALLOC(self->records, sizeof(fasta_record_t) * *cap);
  }
  fasta_record_t *r = self->records + self->n_records++;
  r->name = xstrdup(name);
  r->description = xstrdup(description);
  XMALLOC(r->seq, seq->len + 1);
  memcpy(r->seq, seq->s, seq->len);
  r->seq[seq->len] = '\0';
  r->seq_len = seq->len;
}

/* split a header line (without the leading '>') into a name (the first
   word) and a description (the rest, without leading spaces) */
static void split_header(char *line, char **name, char **description) {
  *name = line;
  char *p = line;
  while (*p != '\0' && !isspace((unsigned char)*p))
    p++;
  if (*p == '\0') {
    *description = p;
    return;
  }
  *p++ = '\0';
  while (*p != '\0' && isspace((unsigned char)*p))
    p++;
  *description = p;
}

fasta_t* fasta_read(const char *path, Logger *logger) {
  assert(path != NULL);
  FILE *file = xfopen(path, NULL, "r", "reference FASTA", logger);
  char *data;
  size_t size;
  read_file(file, &data, &size, logger);
  xfclose(&file, logger);
  fasta_t *self;
  XMALLOC(self, sizeof(fasta_t));
  self->records = NULL;
  self->n_records = 0;
  size_t cap = 0;
  char *name = NULL, *description = NULL;
  bool in_record = false;
  strbuf_t seq;
  strbuf_init(&seq);
  char *line = data, *end = data + size;
  while (line < end) {
    char *nl = memchr(line, '\n', (size_t)(end - line));
    char *line_end = nl != NULL ? nl : end;
    /* strip a trailing carriage return (CRLF files) */
    if (line_end > line && line_end[-1] == '\r')
      line_end[-1] = '\0';
    else
      *line_end = '\0';
    if (line[0] == '>') {
      if (in_record)
        add_record(self, &cap, name, description, &seq, path, logger);
      split_header(line + 1, &name, &description);
      LCHECK(name[0] != '\0', logger,
             "%s: a record has an empty name\n", path);
      strbuf_clear(&seq);
      in_record = true;
    } else if (line[0] != '\0') {
      LCHECK(in_record, logger,
             "%s: sequence before the first header line\n", path);
      for (char *p = line; *p != '\0'; p++) {
        if (isspace((unsigned char)*p))
          continue;
        char c = (char)toupper((unsigned char)*p);
        LCHECK(primers_iupac_bases(c) != NULL, logger,
               "%s: invalid base '%c' in \"%s\"\n", path, *p, name);
        strbuf_append_char(&seq, c);
      }
    }
    line = line_end + 1;
  }
  if (in_record)
    add_record(self, &cap, name, description, &seq, path, logger);
  strbuf_free(&seq);
  free(data);
  LCHECK(self->n_records > 0, logger,
         "%s: no sequence found (not a FASTA file?)\n", path);
  return self;
}

void fasta_delete(fasta_t *self) {
  if (self == NULL)
    return;
  for (size_t i = 0; i < self->n_records; i++) {
    free(self->records[i].name);
    free(self->records[i].description);
    free(self->records[i].seq);
  }
  free(self->records);
  free(self);
}

/* write <content> to a new file <path> */
static void write_file(const char *path, const char *content,
    Logger *logger) {
  FILE *f = xfopen(path, NULL, "w", "test FASTA", logger);
  fputs(content, f);
  xfclose(&f, logger);
}

void fasta_test(const char *tmpdir) {
  Logger *logger = new_logger(stderr, LOG_ERROR);
  char *path = path_join(tmpdir, "test.fasta");
  write_file(path,
      ">seg1 H5 HA segment\nACGTACGT\nACGT\n"
      ">seg2\r\nacgtN\r\n", logger);
  fasta_t *fa = fasta_read(path, logger);
  test_eq_size(fa->n_records, 2);
  test_eq_str(fa->records[0].name, "seg1");
  test_eq_str(fa->records[0].description, "H5 HA segment");
  test_eq_str(fa->records[0].seq, "ACGTACGTACGT");
  test_eq_size(fa->records[0].seq_len, 12);
  /* lower case is upper-cased, CRLF handled, description may be empty */
  test_eq_str(fa->records[1].name, "seg2");
  test_eq_str(fa->records[1].description, "");
  test_eq_str(fa->records[1].seq, "ACGTN");
  fasta_delete(fa);
  xfunlink(path, NULL, logger);
  free(path);
  delete_logger(logger);
}
