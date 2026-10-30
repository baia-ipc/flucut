/* (c) 2026 Giorgio Gonnella */

/*
   Example program of the libcbaia manual: count the reads and the bases
   of FASTQ or BAM files (or of the files of reads in directories).
*/

#define VERSION "1.0"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "alloc.h"
#include "argparse-def.h"
#include "file.h"
#include "logger.h"
#include "reads-reader.h"

#define DEFAULT_MIN_LENGTH 0
#define BATCH_SIZE 1000

// %-PURPOSE-%
#define PURPOSE "Count the reads and bases of FASTQ or BAM files\n"\
  "\n"\
  "Usage:\n"\
  "   %s [options] <input>...\n"\
  "\n"\
  "Arguments:\n"\
  "   <input>  FASTQ (optionally gzipped) or BAM file, or a directory\n"\
  "            (its files of reads are counted)\n"\
  "\n"\
  "Options:\n"\
  "   -m, --min-length <n>  count only the reads of at least <n> bases\n"\
  "                         (default: " XSTR(DEFAULT_MIN_LENGTH) ")\n"\
  "\n"\
  COMMON_OPTS_HELP

void parseargs(int argc, char *argv[],
    char ***inputs,
    size_t *n_inputs,
    long *min_length,
    Logger **logger)
{
  PARSE_START;
  int c = 0;
  *min_length = DEFAULT_MIN_LENGTH;
  static struct option longopts[] = {
    {"min-length", required_argument, NULL, 'm'},
    COMMON_LONGOPTS
    LONGOPTS_END
  };
  opterr = 0;
  while ((c = getopt_long(argc, argv, COMMON_OPTS_STR "m:", longopts,
                          NULL)) != -1) {
    switch (c) {
      case 'm':
        PARSE_INTEGER_OPT_MIN('m', min_length, 0);
        break;
      COMMON_OPTS_PARSE_CASES;
    }
  }
  CHECK_MIN_NARGS(1);
  *inputs = argv + optind;
  *n_inputs = (size_t)(argc - optind);
  INIT_LOGGER(*logger);
  PARSE_LOG("Inputs: %zu\n", *n_inputs);
  PARSE_LOG("Minimum read length: %ld\n", *min_length);
}

/* list_dir() filter: the files of reads */
static bool is_reads_file(const char *name, const char *path, void *data) {
  (void)data;
  return reads_reader_is_reads_filename(name) && !is_directory(path);
}

/* count the reads of a file and output a line with the counts */
static void count_file(const char *path, long min_length,
    fastq_batch_t *batch, Logger *logger) {
  reads_reader_t *reader = reads_reader_new(path, logger);
  size_t reads = 0, bases = 0, n;
  while ((n = reads_reader_read_batch(reader, batch, BATCH_SIZE)) > 0)
    for (size_t i = 0; i < n; i++) {
      size_t len = batch->records[i].seq_len;
      if ((long)len >= min_length) {
        reads++;
        bases += len;
      }
    }
  printf("%s\t%s\t%zu\t%zu\n", path,
         reads_reader_is_bam(reader) ? "BAM" : "FASTQ", reads, bases);
  lprintf(logger, LOG_INFO, "%s: %zu reads\n", path, reads);
  reads_reader_delete(reader);
}

int main(int argc, char *argv[]) {
  char **inputs;
  size_t n_inputs;
  long min_length;
  Logger *logger;
  parseargs(argc, argv, &inputs, &n_inputs, &min_length, &logger);
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  printf("file\tformat\treads\tbases\n");
  for (size_t i = 0; i < n_inputs; i++) {
    if (!is_directory(inputs[i])) {
      count_file(inputs[i], min_length, &batch, logger);
      continue;
    }
    size_t n_files;
    char **files = list_dir(inputs[i], is_reads_file, NULL, &n_files,
                            logger);
    if (n_files == 0)
      lprintf(logger, LOG_WARN, "No files of reads in %s\n", inputs[i]);
    for (size_t j = 0; j < n_files; j++) {
      count_file(files[j], min_length, &batch, logger);
      free(files[j]);
    }
    free(files);
  }
  fastq_batch_free(&batch);
  lprintf(logger, LOG_SUCCESS, "%zu inputs processed\n", n_inputs);
  delete_logger(logger);
  return EXIT_SUCCESS;
}
