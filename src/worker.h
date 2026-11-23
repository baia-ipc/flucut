/* (c) 2026 Giorgio Gonnella */

/*
   Trimming of the reads of a file (FASTQ or BAM) by a number of worker
   threads: each thread reads a batch of reads in turn
   (under a lock), trims it with its own objects (trimmer, compressor,
   buffers, counts), and writes its compressed output in the order of
   the input, waiting for its turn, so that the output does not depend on
   the number of threads. The batches have WORKER_BATCH_SIZE reads, each
   compressed (with libdeflate) as a gzip member of its own; the
   concatenation of the members is a valid gzip file, also when a file is
   appended to (live mode). With an output file for each category of
   reads (read-outputs.h), the reads of a batch are compressed apart for
   each file: the order of the reads of each file is that of the input.
*/

#ifndef WORKER_H
#define WORKER_H

#include <stdio.h>
#include <stddef.h>
#include "counter.h"
#include "fastq-reader.h"
#include "libdeflate.h"
#include "logger.h"
#include "read-outputs.h"
#include "strbuf.h"
#include "trimmer.h"

/* reads of a batch */
#define WORKER_BATCH_SIZE 1000
/* compression level of the output: the range of libdeflate, and the
   default (see the --compress-level option) */
#define WORKER_MIN_COMPRESS_LEVEL 1
#define WORKER_MAX_COMPRESS_LEVEL 12
#define WORKER_DEFAULT_COMPRESS_LEVEL 6

/* state shared by the workers during the processing of a file */
typedef struct worker_job_t worker_job_t;

/* state of a worker */
typedef struct {
  worker_job_t *job;
  trimmer_t *trimmer;
  struct libdeflate_compressor *compressor;
  counter_t *counts;
  fastq_batch_t batch;
  /* the trimmed reads of a batch and their compression, for each
     output (see trimmer_n_outputs()) */
  size_t n_outputs;
  strbuf_t text[TRIMMER_MAX_OUTPUTS], gz[TRIMMER_MAX_OUTPUTS];
} worker_t;

/* Initialize a worker, with a trimmer for <config> (which must exist as
   long as the worker is used) and a compressor of the output with
   <compress_level> (WORKER_MIN_COMPRESS_LEVEL to
   WORKER_MAX_COMPRESS_LEVEL) */
void worker_init(worker_t *w, const trimmer_config_t *config,
    int compress_level);

/* Classify the reads only against the references <allowed> (NULL: all);
   set before processing a group (see reference_set_restrict) */
void worker_restrict_references(worker_t *w, const bool *allowed);

/* Free the memory of a worker */
void worker_free(worker_t *w);

/*
   Trim the reads of <reads_file> with the <n_workers> <workers> (in
   threads, if more than one); the trimmed reads are appended to the open
   files of <out> (gzipped FASTQ; as many as the outputs of the trimmers
   of the workers), and the counts are added to <counts>
*/
void worker_process_file(worker_t *workers, size_t n_workers,
    const char *reads_file, read_outputs_t *out, counter_t *counts,
    Logger *logger);

#endif
