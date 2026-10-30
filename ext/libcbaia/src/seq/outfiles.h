/* (c) 2026 Giorgio Gonnella */

/*
   Output files of a program writing the reads to several destinations,
   e.g. a demultiplexer (a destination for each sample), with one
   compressed file for each destination:
   - for an input file: <outdir>/<input basename>.<destination>.fastq.gz,
     where the suffix of the input file (.fastq, .fq, optionally followed
     by .gz, or .bam; see reads_reader_is_reads_filename) is removed;
   - for a group of input files (directory mode):
     <outdir>/<group>/<destination>.fastq.gz.

   The files are created when the first data are written to them; the data
   are written already compressed (gzip members).

   BAM output (see outfiles_set_bam): the files are <...>.bam instead of
   <...>.fastq.gz, the data are BGZF blocks, the header is written to each
   new file and the BGZF end-of-file marker when it is closed (it is removed
   before appending to an existing file).
*/

#ifndef OUTFILES_H
#define OUTFILES_H

#include <stdbool.h>
#include <stddef.h>
#include "logger.h"

#define OUTFILES_COMPRESSLEVEL 6

typedef struct outfiles_t outfiles_t;

/*
   Create the output files handler for the input file <reads_file>; the
   output directory is created if it does not exist
*/
outfiles_t* outfiles_new(const char *outdir, const char *reads_file,
    size_t n_destinations, Logger *logger);

/*
   Create the output files handler for a group of input files; the output
   directory of the group is created if it does not exist
*/
outfiles_t* outfiles_new_group(const char *outdir, const char *group,
    size_t n_destinations, Logger *logger);

/*
   Append the data to existing output files, instead of overwriting them
   (to be set before writing)
*/
void outfiles_set_append(outfiles_t *self, bool append);

/*
   BAM output: the header of the file of a destination (BGZF-compressed;
   copied); to be set for each destination, before writing
*/
void outfiles_set_bam(outfiles_t *self, size_t destination,
    const char *header, size_t len);

/* Suffix of the output files */
#define OUTFILES_FASTQ_SFX ".fastq.gz"
#define OUTFILES_BAM_SFX ".bam"

/* Write compressed data to the file of a destination (with its name) */
void outfiles_write(outfiles_t *self, size_t destination, const char *name,
    const char *data, size_t len);

/* Close all files and delete the handler (NULL is accepted) */
void outfiles_delete(outfiles_t *self);

/* Unit tests */
void outfiles_test(const char *tmpdir);

#endif
