/* (c) 2026 Giorgio Gonnella */

/*
   Reader of the reads of a FASTQ file (optionally gzipped) or of an
   unaligned BAM file, recognized from the content, with the same
   interface (batches of FASTQ records, see fastq-reader.h); functions for
   recognizing the files of reads by their name, counting the reads, and
   checking whether a file is complete (e.g. when it is written by another
   program, such as a basecaller).
*/

#ifndef READS_READER_H
#define READS_READER_H

#include <stdbool.h>
#include <stddef.h>
#include "fastq-reader.h"
#include "logger.h"

typedef struct reads_reader_t reads_reader_t;

/*
   Open a FASTQ or BAM file; the program is exited with an error message
   if the file cannot be read
*/
reads_reader_t* reads_reader_new(const char *filename, Logger *logger);

/* Whether the file read is a BAM file */
bool reads_reader_is_bam(const reads_reader_t *self);

/*
   Read the next (at most) max_records records into the batch (as
   fastq_reader_read_batch()); returns their number (0: end of file)
*/
size_t reads_reader_read_batch(reads_reader_t *self, fastq_batch_t *batch,
    size_t max_records);

/* Close the file and delete the reader (NULL is accepted) */
void reads_reader_delete(reads_reader_t *self);

/*
   Whether a file name is that of a file of reads: suffix .fastq, .fq
   (optionally followed by .gz) or .bam
*/
bool reads_reader_is_reads_filename(const char *name);

/*
   Length of the suffix of a file of reads (see
   reads_reader_is_reads_filename) at the end of the name, 0 if none
*/
size_t reads_reader_filename_suffix_len(const char *name);

/*
   Number of reads of a FASTQ or BAM file (the file is read and parsed, see
   fastq_reader_count() and bam_reader_count()); the program is exited with
   an error message if the file cannot be read or is not valid
*/
size_t reads_reader_count(const char *filename, Logger *logger);

/*
   Whether a file of reads is complete: a gzipped file is not truncated
   and its content ends with a newline; a FASTQ file ends with a newline; a
   BAM file (suffix .bam) ends with the BGZF end-of-file marker; false if
   the file cannot be read. A gzipped file is decompressed completely
   (as its last character is not known otherwise): it takes the time of
   reading the file; an uncompressed or BAM file is read only at its
   start and end.
*/
bool reads_reader_file_complete(const char *path);

/* Unit tests */
void reads_reader_test(const char *tmpdir);

#endif
