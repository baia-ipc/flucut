/* (c) 2026 Giorgio Gonnella */

/*
   The output files of a group of input files (see input-groups.h), to
   which the trimmed reads are written, gzipped:
   - one file, <outdir>/<group>.fastq.gz; or
   - one for each category of reads (read-cover.h; the --by-category
     option): <outdir>/<group>.<category>.fastq.gz, where <category> is
     full, internal, deletion or other. All of them are created, also
     those to which no read is written.
   A file without reads is a valid gzip file (an empty gzip member).
*/

#ifndef READ_OUTPUTS_H
#define READ_OUTPUTS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include "logger.h"
#include "read-cover.h"

typedef struct {
  size_t n;                         /* files: 1 or READ_N_CATEGORIES */
  char *paths[READ_N_CATEGORIES];
  FILE *files[READ_N_CATEGORIES];   /* NULL: not open */
  /* how the output is named in the messages and the tables: the path of
     the file, or <outdir>/<group>.{full,...}.fastq.gz */
  char *label;
} read_outputs_t;

/* The outputs of the group <group> in <outdir>; <by_category>: a file
   for each category of reads */
read_outputs_t* read_outputs_new(const char *outdir, const char *group,
    bool by_category);

/* Open the files, emptied or, if <append>, appended to */
void read_outputs_open(read_outputs_t *self, bool append, Logger *logger);

/* Write <len> bytes of <data> to the file <k> (with by category, the
   category of the reads) */
void read_outputs_write(read_outputs_t *self, size_t k, const char *data,
    size_t len, Logger *logger);

/* Close the files; an empty gzip member is written to those still
   without any byte */
void read_outputs_close(read_outputs_t *self, Logger *logger);

/* Remove the files, if they exist */
void read_outputs_remove(const read_outputs_t *self, Logger *logger);

/* Whether <path> is the path of one of the files */
bool read_outputs_has_path(const read_outputs_t *self, const char *path);

/* Delete (NULL is accepted) */
void read_outputs_delete(read_outputs_t *self);

/* Unit tests (<tmpdir>: directory for temporary files, left empty) */
void read_outputs_test(const char *tmpdir);

#endif
