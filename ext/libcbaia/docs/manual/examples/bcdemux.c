/* (c) 2026 Giorgio Gonnella */

/*
   Example program of the libcbaia manual: demultiplex the reads of a
   FASTQ or BAM file by a barcode at the start of the reads, corrected
   using a whitelist; the reads of each barcode, without the barcode, are
   written to <outdir>/<input>.<barcode>.fastq.gz, the others to
   <outdir>/<input>.unassigned.fastq.gz.
*/

#define VERSION "1.0"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "alloc.h"
#include "argparse-def.h"
#include "fastq-writer.h"
#include "gzip-compressor.h"
#include "khash.h"
#include "logger.h"
#include "outfiles.h"
#include "reads-reader.h"
#include "strbuf.h"
#include "whitelist.h"

#define DEFAULT_MAX_DIST 1
#define DEFAULT_OUTDIR "."
#define BATCH_SIZE 10000

// %-PURPOSE-%
#define PURPOSE "Demultiplex reads by a barcode at their start\n"\
  "\n"\
  "Usage:\n"\
  "   %s [options] -w <whitelist> <reads>\n"\
  "\n"\
  "Arguments:\n"\
  "   <reads>         FASTQ (optionally gzipped) or BAM file\n"\
  "\n"\
  "Options:\n"\
  "   -w, --whitelist <file>  barcodes, one per line (required)\n"\
  "   -d, --max-dist <n>      maximum edit distance of the barcodes, 0 to 3\n"\
  "                           (default: " XSTR(DEFAULT_MAX_DIST) ")\n"\
  "   -o, --outdir <dir>      output directory (default: " DEFAULT_OUTDIR ")\n"\
  "\n"\
  COMMON_OPTS_HELP

void parseargs(int argc, char *argv[],
    char **reads_filename,
    char **whitelist_filename,
    long *max_dist,
    char **outdir,
    Logger **logger)
{
  PARSE_START;
  int c = 0;
  *whitelist_filename = NULL;
  *max_dist = DEFAULT_MAX_DIST;
  *outdir = DEFAULT_OUTDIR;
  static struct option longopts[] = {
    {"whitelist", required_argument, NULL, 'w'},
    {"max-dist", required_argument, NULL, 'd'},
    {"outdir", required_argument, NULL, 'o'},
    COMMON_LONGOPTS
    LONGOPTS_END
  };
  opterr = 0;
  while ((c = getopt_long(argc, argv, COMMON_OPTS_STR "w:d:o:", longopts,
                          NULL)) != -1) {
    switch (c) {
      case 'w':
        *whitelist_filename = optarg;
        break;
      case 'd':
        PARSE_INTEGER_OPT_MIN_MAX('d', max_dist, 0, WHITELIST_MAX_DIST);
        break;
      case 'o':
        *outdir = optarg;
        break;
      COMMON_OPTS_PARSE_CASES;
    }
  }
  if (*whitelist_filename == NULL)
    PARSE_ERR("the whitelist (-w) is required");
  CHECK_NARGS(1);
  *reads_filename = argv[optind];
  INIT_LOGGER(*logger);
  PARSE_LOG("Reads: %s\n", *reads_filename);
  PARSE_LOG("Whitelist: %s\n", *whitelist_filename);
  PARSE_LOG("Maximum barcode distance: %ld\n", *max_dist);
  PARSE_LOG("Output directory: %s\n", *outdir);
}

/* encoded barcode -> number of the destination */
KHASH_MAP_INIT_INT64(dest, size_t)

int main(int argc, char *argv[]) {
  char *reads_filename, *whitelist_filename, *outdir;
  long max_dist;
  Logger *logger;
  parseargs(argc, argv, &reads_filename, &whitelist_filename, &max_dist,
            &outdir, &logger);

  /* the destinations: the barcodes of the whitelist, then "unassigned" */
  whitelist_t *wl = whitelist_new(whitelist_filename, logger);
  uint64_t *codes;
  size_t n_barcodes = whitelist_codes(wl, &codes);
  size_t blen = whitelist_barcode_length(wl), unassigned = n_barcodes;
  char **names;
  XMALLOC(names, sizeof(char *) * (n_barcodes + 1));
  /* khash does not exit on allocation failures: kh_init() returns NULL,
     kh_put() sets absent to -1 */
  khash_t(dest) *dest = kh_init(dest);
  LCHECK(dest != NULL, logger, "Failed allocating the destinations\n");
  for (size_t i = 0; i < n_barcodes; i++) {
    XMALLOC(names[i], blen + 1);
    whitelist_decode(wl, codes[i], names[i]);
    int absent;
    khiter_t k = kh_put(dest, dest, codes[i], &absent);
    LCHECK(absent != -1, logger, "Failed allocating the destinations\n");
    kh_value(dest, k) = i;
  }
  names[unassigned] = "unassigned";

  outfiles_t *outfiles = outfiles_new(outdir, reads_filename,
                                      n_barcodes + 1, logger);
  gzip_compressor_t *compressor =
    gzip_compressor_new(OUTFILES_COMPRESSLEVEL);
  strbuf_t *fastq, compressed;
  XMALLOC(fastq, sizeof(strbuf_t) * (n_barcodes + 1));
  for (size_t d = 0; d <= n_barcodes; d++)
    strbuf_init(&fastq[d]);
  strbuf_init(&compressed);
  size_t *counts;
  XCALLOC(counts, n_barcodes + 1, sizeof(size_t));

  reads_reader_t *reader = reads_reader_new(reads_filename, logger);
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  size_t n;
  while ((n = reads_reader_read_batch(reader, &batch, BATCH_SIZE)) > 0) {
    /* the records of each destination, in FASTQ format */
    for (size_t i = 0; i < n; i++) {
      const fastq_record_t *r = &batch.records[i];
      const char *data = batch.data.s;
      whitelist_match_t m;
      whitelist_match(wl, data + r->seq_from, r->seq_len, false, max_dist,
                      &m);
      size_t d = unassigned, skip = 0;
      if (m.status == WHITELIST_EXACT || m.status == WHITELIST_CORRECTED) {
        d = kh_value(dest, kh_get(dest, dest, m.barcode));
        skip = m.consumed;
      }
      fastq_writer_record(&fastq[d], data + r->title_from, r->title_len,
                          data + r->seq_from + skip,
                          data + r->qual_from + skip, r->seq_len - skip);
      counts[d]++;
    }
    /* compressed, and written to the files of the destinations */
    for (size_t d = 0; d <= n_barcodes; d++) {
      if (fastq[d].len == 0)
        continue;
      strbuf_clear(&compressed);
      gzip_compressor_run(compressor, fastq[d].s, fastq[d].len,
                          &compressed);
      outfiles_write(outfiles, d, names[d], compressed.s, compressed.len);
      strbuf_clear(&fastq[d]);
    }
  }

  for (size_t d = 0; d <= n_barcodes; d++)
    if (counts[d] > 0)
      lprintf(logger, LOG_SUCCESS, "%s: %zu reads\n", names[d], counts[d]);

  fastq_batch_free(&batch);
  reads_reader_delete(reader);
  outfiles_delete(outfiles);
  gzip_compressor_delete(compressor);
  for (size_t d = 0; d <= n_barcodes; d++)
    strbuf_free(&fastq[d]);
  strbuf_free(&compressed);
  free(fastq);
  free(counts);
  for (size_t i = 0; i < n_barcodes; i++)
    free(names[i]);
  free(names);
  free(codes);
  kh_destroy(dest, dest);
  whitelist_delete(wl);
  delete_logger(logger);
  return EXIT_SUCCESS;
}
