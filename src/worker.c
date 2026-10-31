/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "alloc.h"
#include "reads-reader.h"
#include "worker.h"

/* state shared by the threads */
struct worker_job_t {
  reads_reader_t *reader;
  FILE *out;
  const char *out_path;
  pthread_mutex_t input_mutex, output_mutex;
  pthread_cond_t output_cond;
  size_t next_batch, next_to_write;
  Logger *logger;
};

void worker_init(worker_t *w, const trimmer_config_t *config,
    int compress_level) {
  assert(w != NULL);
  assert(config != NULL);
  assert(compress_level >= WORKER_MIN_COMPRESS_LEVEL &&
         compress_level <= WORKER_MAX_COMPRESS_LEVEL);
  w->job = NULL;
  w->trimmer = trimmer_new(config);
  w->compressor = libdeflate_alloc_compressor(compress_level);
  if (w->compressor == NULL)
    fatal_printf("Initialization of the gzip compression failed\n");
  w->counts = counter_new(config->scheme, config->references);
  fastq_batch_init(&w->batch);
  strbuf_init(&w->text);
  strbuf_init(&w->gz);
}

void worker_restrict_subtype(worker_t *w, int subtype) {
  trimmer_restrict_subtype(w->trimmer, subtype);
}

void worker_free(worker_t *w) {
  trimmer_delete(w->trimmer);
  libdeflate_free_compressor(w->compressor);
  counter_delete(w->counts);
  fastq_batch_free(&w->batch);
  strbuf_free(&w->text);
  strbuf_free(&w->gz);
}

/*
   Compress <len> bytes of <data> as a gzip member, appended to <out>
   (in a space of the size of the bound of libdeflate: the compression
   cannot fail for lack of space)
*/
static void worker_compress(struct libdeflate_compressor *compressor,
    const char *data, size_t len, strbuf_t *out) {
  size_t bound = libdeflate_gzip_compress_bound(compressor, len);
  char *dest = strbuf_reserve(out, bound);
  size_t n = libdeflate_gzip_compress(compressor, data, len, dest, bound);
  if (n == 0)
    fatal_printf("gzip compression failed\n");
  strbuf_commit(out, n);
}

/*
   Trim batches of reads until the end of the input: the batches are read
   in turn, trimmed in parallel, and their output is written in the input
   order; the counts are kept by each worker, and summed by the caller
*/
static void *worker_run(void *data) {
  worker_t *w = data;
  worker_job_t *job = w->job;
  while (true) {
    pthread_mutex_lock(&job->input_mutex);
    size_t batch_num = job->next_batch++;
    size_t n = reads_reader_read_batch(job->reader, &w->batch,
                                       WORKER_BATCH_SIZE);
    pthread_mutex_unlock(&job->input_mutex);
    if (n == 0)
      break;
    strbuf_clear(&w->text);
    strbuf_clear(&w->gz);
    trimmer_run(w->trimmer, &w->batch, 0, n, &w->text, w->counts);
    if (w->text.len > 0)
      worker_compress(w->compressor, w->text.s, w->text.len, &w->gz);
    pthread_mutex_lock(&job->output_mutex);
    while (job->next_to_write != batch_num)
      pthread_cond_wait(&job->output_cond, &job->output_mutex);
    LCHECK(fwrite(w->gz.s, 1, w->gz.len, job->out) == w->gz.len,
           job->logger, "%s: write error\n", job->out_path);
    job->next_to_write++;
    pthread_cond_broadcast(&job->output_cond);
    pthread_mutex_unlock(&job->output_mutex);
  }
  return NULL;
}

void worker_process_file(worker_t *workers, size_t n_workers,
    const char *reads_file, FILE *out, const char *out_path,
    counter_t *counts, Logger *logger) {
  assert(workers != NULL);
  assert(n_workers > 0);
  assert(counts != NULL);
  worker_job_t job;
  job.reader = reads_reader_new(reads_file, logger);
  job.out = out;
  job.out_path = out_path;
  job.next_batch = 0;
  job.next_to_write = 0;
  job.logger = logger;
  pthread_mutex_init(&job.input_mutex, NULL);
  pthread_mutex_init(&job.output_mutex, NULL);
  pthread_cond_init(&job.output_cond, NULL);
  for (size_t t = 0; t < n_workers; t++) {
    workers[t].job = &job;
    counter_clear(workers[t].counts);
  }
  if (n_workers == 1)
    worker_run(workers);
  else {
    pthread_t *threads;
    XMALLOC(threads, sizeof(pthread_t) * n_workers);
    for (size_t t = 0; t < n_workers; t++)
      LCHECK(pthread_create(threads + t, NULL, worker_run, workers + t)
             == 0, logger, "Cannot start a thread\n");
    for (size_t t = 0; t < n_workers; t++)
      pthread_join(threads[t], NULL);
    free(threads);
  }
  /* (the counts are sums: their order does not matter) */
  for (size_t t = 0; t < n_workers; t++) {
    counter_merge(counts, workers[t].counts);
    workers[t].job = NULL;
  }
  pthread_mutex_destroy(&job.input_mutex);
  pthread_mutex_destroy(&job.output_mutex);
  pthread_cond_destroy(&job.output_cond);
  reads_reader_delete(job.reader);
}
