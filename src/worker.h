/* (c) 2026 Giorgio Gonnella */

/*
   Trimming of the reads of a file (FASTQ or BAM) by a number of worker
   threads (from Lumipore): each thread reads a batch of reads in turn
   (under a lock), trims it with its own objects (trimmer, compressor,
   buffers, counts), and writes its compressed output in the order of
   the input, waiting for its turn, so that the output does not depend on
   the number of threads. The batches have WORKER_BATCH_SIZE reads, each
   compressed as a gzip member of its own.
*/

#ifndef WORKER_H
#define WORKER_H

#include <stdio.h>
#include <stddef.h>
#include "counter.h"
#include "fastq-reader.h"
#include "gzip-compressor.h"
#include "logger.h"
#include "strbuf.h"
#include "trimmer.h"

/* reads of a batch */
#define WORKER_BATCH_SIZE 1000
/* compression level of the output */
#define WORKER_COMPRESS_LEVEL 6

/* state shared by the workers during the processing of a file */
typedef struct worker_job_t worker_job_t;

/* state of a worker */
typedef struct {
  worker_job_t *job;
  trimmer_t *trimmer;
  gzip_compressor_t *compressor;
  counter_t *counts;
  fastq_batch_t batch;
  strbuf_t text, gz;
} worker_t;

/* Initialize a worker, with a trimmer for <config> (which must exist as
   long as the worker is used) */
void worker_init(worker_t *w, const trimmer_config_t *config);

/* Free the memory of a worker */
void worker_free(worker_t *w);

/*
   Trim the reads of <reads_file> with the <n_workers> <workers> (in
   threads, if more than one); the trimmed reads are appended to <out>
   (gzipped FASTQ, path <out_path>, for the messages), and the counts are
   added to <counts>
*/
void worker_process_file(worker_t *workers, size_t n_workers,
    const char *reads_file, FILE *out, const char *out_path,
    counter_t *counts, Logger *logger);

#endif
