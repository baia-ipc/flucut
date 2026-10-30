/* (c) 2026 Giorgio Gonnella */

/*
   Reading of FASTQ files (optionally gzip-compressed) in batches of
   records.

   The parsing follows the Biopython FastqGeneralIterator: the title line is
   kept as it is (except trailing whitespace); sequences and qualities can
   span several lines (a line starting with '@' ends the qualities, if they
   are at least as long as the sequence); the title on the "+" line is
   optional, but must be equal to the title of the record, if present.
   Unlike Biopython, which skips any text before the first record, a file
   whose first line does not start with '@' is not valid.
*/

#ifndef FASTQ_READER_H
#define FASTQ_READER_H

#include <stdbool.h>
#include <stddef.h>
#include "logger.h"
#include "strbuf.h"

typedef struct fastq_reader_t fastq_reader_t;

/* a record: offsets and lengths in the data of the batch */
typedef struct {
  size_t title_from, title_len;
  size_t seq_from, seq_len;
  size_t qual_from;          /* the qualities have length seq_len */
  bool qual_missing;         /* no qualities (BAM input): given as '!' */
} fastq_record_t;

/* a batch of records */
typedef struct {
  strbuf_t data;
  fastq_record_t *records;
  size_t n_records, records_cap;
} fastq_batch_t;

/* Open a FASTQ file; exit with an error message on failure */
fastq_reader_t* fastq_reader_new(const char *filename, Logger *logger);

/*
   Read up to <max_records> records into <batch> (its previous content is
   removed); returns the number of records read (0 at the end of the file).
   The program is exited with an error message if the file is not valid.
*/
size_t fastq_reader_read_batch(fastq_reader_t *self, fastq_batch_t *batch,
    size_t max_records);

/* Close the file and delete the reader (NULL is accepted) */
void fastq_reader_delete(fastq_reader_t *self);

/*
   Number of records of a FASTQ file (the file is read, as by
   fastq_reader_read_batch()); the program is exited with an error message
   if the file cannot be read or is not valid
*/
size_t fastq_reader_count(const char *filename, Logger *logger);

/* Initialize an empty batch */
void fastq_batch_init(fastq_batch_t *batch);

/* Free the memory of a batch */
void fastq_batch_free(fastq_batch_t *batch);

/* Unit tests (<tmpdir>: directory for temporary files) */
void fastq_reader_test(const char *tmpdir);

#endif
