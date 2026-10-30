/* (c) 2026 Giorgio Gonnella */

/*
   Reading of BAM files, e.g. the unaligned BAM files written by Dorado
   and MinKNOW, in batches of FASTQ-like records (see fastq-reader.h).

   BAM is a series of BGZF blocks, which are gzip members: the file is
   decompressed with zlib, without htslib.

   For each record:
   - the title is the read name, followed by the tags (SAM text format,
     separated by tabs, as in the FASTQ output of Dorado), except array
     tags (type B, e.g. the move table) and the base modification tags
     (MM, ML, MN), whose positions refer to the whole read;
   - the sequence and qualities are in the orientation of the original
     read (reverse complemented, if the record is mapped to the minus
     strand); quality values above 93 are written as '~' (the maximum of
     SAM); missing qualities are written as '!', and the record is
     marked (qual_missing), so that, written again (ubam-writer.h, with
     NULL qualities), they are missing qualities again;
   - secondary and supplementary alignments are skipped.
*/

#ifndef BAM_READER_H
#define BAM_READER_H

#include <stdbool.h>
#include <stddef.h>
#include "fastq-reader.h"
#include "logger.h"

typedef struct bam_reader_t bam_reader_t;

/* Is the file a BAM file (by its content); exit with an error message if
   it cannot be opened */
bool bam_reader_is_bam(const char *filename, Logger *logger);

/* Open a BAM file and read its header; exit with an error message on
   failure */
bam_reader_t* bam_reader_new(const char *filename, Logger *logger);

/* SAM header text of the file (e.g. its @RG and @PG lines), valid as long
   as the reader exists */
const char *bam_reader_header(const bam_reader_t *self);

/*
   Read up to <max_records> records into <batch> (its previous content is
   removed); returns the number of records read (0 at the end of the file).
   The program is exited with an error message if the file is not valid.
*/
size_t bam_reader_read_batch(bam_reader_t *self, fastq_batch_t *batch,
    size_t max_records);

/* Close the file and delete the reader (NULL is accepted) */
void bam_reader_delete(bam_reader_t *self);

/* Does the file end with the BGZF end-of-file marker, written when a BAM
   file is closed (i.e. the file is complete) */
bool bam_reader_has_eof_marker(const char *filename);

/*
   Number of records of a BAM file (as returned by the reader); the program
   is exited with an error message if the file cannot be read or is not
   valid
*/
size_t bam_reader_count(const char *filename, Logger *logger);

/* Unit tests (<tmpdir>: directory for temporary files) */
void bam_reader_test(const char *tmpdir);

#endif
