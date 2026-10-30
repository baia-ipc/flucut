/* (c) 2026 Giorgio Gonnella */

/*
   Reading of a FASTA file into memory.

   The reference sequences of a run (the segments against which the reads
   are aligned, see reference.h) are a small FASTA file, read whole. Each
   record is a header line starting with '>' (the name is the first word,
   the rest of the line is the description) followed by one or more
   sequence lines. The sequence is stored without newlines and in upper
   case; it may contain IUPAC codes.

   Only plain (uncompressed) FASTA is read: the references are small and
   are not compressed. The reads of a run are read with the FASTQ/BAM
   reader of libcbaia, not with this module.
*/

#ifndef FASTA_H
#define FASTA_H

#include <stddef.h>
#include "logger.h"

typedef struct {
  char *name;          /* first word of the header (without '>') */
  char *description;   /* rest of the header line (may be empty) */
  char *seq;           /* sequence, upper case, 0-terminated */
  size_t seq_len;
} fasta_record_t;

typedef struct {
  fasta_record_t *records;
  size_t n_records;
} fasta_t;

/*
   Read a FASTA file; the program is exited with an error message if the
   file cannot be read or is not valid (no record, a sequence with an
   invalid character, an empty sequence).
*/
fasta_t* fasta_read(const char *path, Logger *logger);

/* Delete (NULL is accepted) */
void fasta_delete(fasta_t *self);

/* Unit tests (<tmpdir>: directory for temporary files) */
void fasta_test(const char *tmpdir);

#endif
