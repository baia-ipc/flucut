/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "alloc.h"
#include "bam-reader.h"
#include "file.h"
#include "sam-header.h"
#include "test.h"
#include "test-bam.h"

void sam_header_inputs_init(sam_header_inputs_t *self) {
  assert(self != NULL);
  strbuf_init(&self->rg_lines);
  strbuf_init(&self->pg_lines);
  self->rg_ids = NULL;
  self->n_rg = 0;
  self->last_pg = NULL;
}

void sam_header_inputs_clear(sam_header_inputs_t *self) {
  assert(self != NULL);
  strbuf_clear(&self->rg_lines);
  strbuf_clear(&self->pg_lines);
  for (size_t i = 0; i < self->n_rg; i++)
    free(self->rg_ids[i]);
  self->n_rg = 0;
  free(self->last_pg);
  self->last_pg = NULL;
}

void sam_header_inputs_free(sam_header_inputs_t *self) {
  if (self == NULL)
    return;
  sam_header_inputs_clear(self);
  free(self->rg_ids);
  self->rg_ids = NULL;
  strbuf_free(&self->rg_lines);
  strbuf_free(&self->pg_lines);
}

const char *sam_header_line_id(const char *line, size_t line_len,
    size_t *len) {
  assert(line != NULL || line_len == 0);
  assert(len != NULL);
  for (const char *p = line; p + 4 <= line + line_len; p++)
    if (p[0] == '\t' && strncmp(p + 1, "ID:", 3) == 0) {
      const char *id = p + 4, *end = id;
      while (end < line + line_len && *end != '\t')
        end++;
      *len = (size_t)(end - id);
      return id;
    }
  return NULL;
}

/* is there a line of the same type (e.g. "@PG") with the same ID in
   the header lines <lines> */
static bool has_id(const strbuf_t *lines, const char *type, const char *id,
    size_t id_len) {
  if (lines->len == 0)  /* (the buffer can be NULL) */
    return false;
  const char *p = lines->s, *end = lines->s + lines->len;
  while (p < end) {
    const char *eol = memchr(p, '\n', (size_t)(end - p));
    if (eol == NULL)
      eol = end;
    size_t len;
    const char *other = sam_header_line_id(p, (size_t)(eol - p), &len);
    if (strncmp(p, type, 3) == 0 && other != NULL && len == id_len &&
        memcmp(other, id, len) == 0)
      return true;
    p = eol + 1;
  }
  return false;
}

void sam_header_inputs_add(sam_header_inputs_t *self, const char *text) {
  assert(self != NULL);
  assert(text != NULL);
  const char *p = text, *end = text + strlen(text);
  while (p < end) {
    const char *eol = memchr(p, '\n', (size_t)(end - p));
    if (eol == NULL)
      eol = end;
    size_t id_len, len = (size_t)(eol - p);
    const char *id = sam_header_line_id(p, len, &id_len);
    bool rg = strncmp(p, "@RG\t", 4) == 0, pg = strncmp(p, "@PG\t", 4) == 0;
    strbuf_t *lines = rg ? &self->rg_lines : &self->pg_lines;
    if (id != NULL && (rg || pg) && !has_id(lines, p, id, id_len)) {
      strbuf_append(lines, p, len);
      strbuf_append_char(lines, '\n');
      char *copy;
      XMALLOC(copy, id_len + 1);
      memcpy(copy, id, id_len);
      copy[id_len] = '\0';
      if (rg) {
        XREALLOC(self->rg_ids, sizeof(char *) * (self->n_rg + 1));
        self->rg_ids[self->n_rg++] = copy;
      } else {
        free(self->last_pg);
        self->last_pg = copy;
      }
    }
    p = eol + 1;
  }
}

void sam_header_inputs_add_files(sam_header_inputs_t *self,
    char *const *files, size_t n_files, Logger *logger) {
  assert(self != NULL);
  assert(files != NULL || n_files == 0);
  for (size_t f = 0; f < n_files; f++) {
    if (!bam_reader_is_bam(files[f], logger))
      continue;
    bam_reader_t *reader = bam_reader_new(files[f], logger);
    sam_header_inputs_add(self, bam_reader_header(reader));
    bam_reader_delete(reader);
  }
}

void sam_header_append_read_groups(strbuf_t *text,
    const sam_header_inputs_t *inputs, const char *name) {
  assert(text != NULL);
  assert(inputs != NULL);
  if (name == NULL) {
    strbuf_append(text, inputs->rg_lines.s, inputs->rg_lines.len);
    return;
  }
  strbuf_append_str(text, "@RG\tID:");
  strbuf_append_str(text, name);
  strbuf_append_str(text, "\tSM:");
  strbuf_append_str(text, name);
  strbuf_append_char(text, '\n');
  if (inputs->rg_lines.len == 0)  /* (the buffer can be NULL) */
    return;
  const char *p = inputs->rg_lines.s,
             *end = inputs->rg_lines.s + inputs->rg_lines.len;
  while (p < end) {
    const char *eol = memchr(p, '\n', (size_t)(end - p));
    strbuf_append_str(text, "@RG");
    /* the fields, with the ID prefixed by the name and SM replaced */
    for (const char *f = p + 3; f < eol; ) {
      const char *fend = f + 1;
      while (fend < eol && *fend != '\t')
        fend++;
      if (strncmp(f, "\tID:", 4) == 0) {
        strbuf_append_str(text, "\tID:");
        strbuf_append_str(text, name);
        strbuf_append_char(text, '.');
        strbuf_append(text, f + 4, (size_t)(fend - f - 4));
      } else if (strncmp(f, "\tSM:", 4) != 0)
        strbuf_append(text, f, (size_t)(fend - f));
      f = fend;
    }
    strbuf_append_str(text, "\tSM:");
    strbuf_append_str(text, name);
    strbuf_append_char(text, '\n');
    p = eol + 1;
  }
}

void sam_header_append_programs(strbuf_t *text,
    const sam_header_inputs_t *inputs, const char *program,
    const char *version, const char *command_line) {
  assert(text != NULL);
  assert(inputs != NULL);
  assert(program != NULL);
  assert(version != NULL);
  assert(command_line != NULL);
  strbuf_append(text, inputs->pg_lines.s, inputs->pg_lines.len);
  /* the ID must be unique (e.g. prog.1 if an input was written by prog) */
  size_t id_size = strlen(program) + 24;
  char *id;
  XMALLOC(id, id_size);
  snprintf(id, id_size, "%s", program);
  for (int i = 1; has_id(&inputs->pg_lines, "@PG", id, strlen(id)); i++)
    snprintf(id, id_size, "%s.%d", program, i);
  strbuf_append_str(text, "@PG\tID:");
  strbuf_append_str(text, id);
  strbuf_append_str(text, "\tPN:");
  strbuf_append_str(text, program);
  if (inputs->last_pg != NULL) {
    strbuf_append_str(text, "\tPP:");
    strbuf_append_str(text, inputs->last_pg);
  }
  strbuf_append_str(text, "\tVN:");
  strbuf_append_str(text, version);
  strbuf_append_str(text, "\tCL:");
  strbuf_append_str(text, command_line);
  strbuf_append_char(text, '\n');
  free(id);
}

char *sam_header_command_line(int argc, char *const argv[]) {
  assert(argv != NULL || argc == 0);
  strbuf_t cl;
  strbuf_init(&cl);
  for (int i = 0; i < argc; i++) {
    if (i > 0)
      strbuf_append_char(&cl, ' ');
    for (const char *p = argv[i]; *p != '\0'; p++)
      strbuf_append_char(&cl, *p == '\t' || *p == '\n' ? ' ' : *p);
  }
  strbuf_append_char(&cl, '\0');
  return cl.s;
}

/* --- unit tests --- */

static void test_text(strbuf_t *text, const char *exp) {
  strbuf_append_char(text, '\0');
  test_eq_str(text->s, exp);
  strbuf_clear(text);
}

/* write a BAM file with the header <text> (no reads) */
static char *write_bam(const char *tmpdir, const char *name,
    const char *text) {
  return test_bam_file(tmpdir, name, text, NULL, 0, true);
}

/* the read groups and programs of the BAM files among the inputs */
static void test_add_files(const char *tmpdir) {
  char *files[3];
  files[0] = write_bam(tmpdir, "cbaia-tests-sam-header1.bam",
                       "@HD\tVN:1.6\n@RG\tID:r1\tSM:x\n"
                       "@PG\tID:basecaller\tPN:dorado\n");
  files[1] = path_join(tmpdir, "cbaia-tests-sam-header.fastq");
  FILE *f = fopen(files[1], "w");
  test_assert(f != NULL);
  fputs("@r1\nACGT\n+\nIIII\n", f);
  fclose(f);
  files[2] = write_bam(tmpdir, "cbaia-tests-sam-header2.bam",
                       "@RG\tID:r1\tSM:x\n@RG\tID:r2\n"
                       "@PG\tID:prog\tPN:prog\tPP:basecaller\n");
  sam_header_inputs_t in;
  sam_header_inputs_init(&in);
  sam_header_inputs_add_files(&in, files, 3, NULL);
  test_eq_size(in.n_rg, 2);
  test_eq_str(in.rg_ids[0], "r1");
  test_eq_str(in.rg_ids[1], "r2");
  test_eq_str(in.last_pg, "prog");
  const char *rg = "@RG\tID:r1\tSM:x\n@RG\tID:r2\n";
  test_assert(in.rg_lines.len == strlen(rg) &&
              memcmp(in.rg_lines.s, rg, strlen(rg)) == 0);
  /* no files */
  sam_header_inputs_clear(&in);
  sam_header_inputs_add_files(&in, NULL, 0, NULL);
  test_eq_size(in.n_rg, 0);
  test_assert(in.last_pg == NULL);
  sam_header_inputs_free(&in);
  for (size_t i = 0; i < 3; i++) {
    unlink(files[i]);
    free(files[i]);
  }
}

void sam_header_test(const char *tmpdir) {
  test_add_files(tmpdir);
  size_t len;
  const char *id = sam_header_line_id("@RG\tID:run1\tSM:x", 16, &len);
  test_assert(id != NULL && len == 4 && strncmp(id, "run1", 4) == 0);
  test_assert(sam_header_line_id("@HD\tVN:1.6", 10, &len) == NULL);
  char *argv[] = {"prog", "-x", "a\tb", "c\nd"};
  char *cl = sam_header_command_line(4, argv);
  test_eq_str(cl, "prog -x a b c d");
  sam_header_inputs_t in;
  sam_header_inputs_init(&in);
  strbuf_t text;
  strbuf_init(&text);
  /* no inputs */
  sam_header_append_read_groups(&text, &in, "S1");
  sam_header_append_programs(&text, &in, "prog", "1.0", cl);
  test_text(&text, "@RG\tID:S1\tSM:S1\n"
            "@PG\tID:prog\tPN:prog\tVN:1.0\tCL:prog -x a b c d\n");
  /* two inputs, with a repeated read group, and an earlier prog */
  sam_header_inputs_add(&in, "@HD\tVN:1.6\n@RG\tID:r1\tSM:x\tPL:ONT\n"
                        "@PG\tID:basecaller\tPN:dorado\n");
  sam_header_inputs_add(&in, "@RG\tID:r1\tSM:x\tPL:ONT\n@RG\tID:r2\n"
                        "@PG\tID:prog\tPN:prog\tPP:basecaller\n");
  test_eq_size(in.n_rg, 2);
  test_eq_str(in.rg_ids[1], "r2");
  test_eq_str(in.last_pg, "prog");
  sam_header_append_read_groups(&text, &in, NULL);
  test_text(&text, "@RG\tID:r1\tSM:x\tPL:ONT\n@RG\tID:r2\n");
  sam_header_append_read_groups(&text, &in, "S1");
  sam_header_append_programs(&text, &in, "prog", "1.0", "prog");
  test_text(&text, "@RG\tID:S1\tSM:S1\n@RG\tID:S1.r1\tPL:ONT\tSM:S1\n"
            "@RG\tID:S1.r2\tSM:S1\n@PG\tID:basecaller\tPN:dorado\n"
            "@PG\tID:prog\tPN:prog\tPP:basecaller\n"
            "@PG\tID:prog.1\tPN:prog\tPP:prog\tVN:1.0\tCL:prog\n");
  sam_header_inputs_clear(&in);
  test_eq_size(in.n_rg, 0);
  test_assert(in.last_pg == NULL);
  strbuf_free(&text);
  sam_header_inputs_free(&in);
  free(cl);
}
