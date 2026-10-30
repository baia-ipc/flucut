/* (c) 2026 Giorgio Gonnella */

/*
   The options of flucut, their help text and their parsing.
*/

#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdbool.h>
#include <stddef.h>
#include "amplicon.h"
#include "logger.h"

#define PROGRAM_NAME "flucut"
#define VERSION "0.1"

/* defaults of the options which are not defaults of a module */
#define OPTIONS_DEFAULT_OUTDIR "flucut-output"
#define OPTIONS_DEFAULT_MAX_ERRORS 0.2
#define OPTIONS_DEFAULT_MIN_LENGTH 1
#define OPTIONS_MAX_THREADS 256

typedef struct {
  /* the inputs (files of reads, or directories of files of reads): the
     arguments of the command line */
  char **inputs;
  size_t n_inputs;
  const char *outdir;
  const char *scheme;
  trim_policy_t policy;
  double max_errors;
  long window;
  const char *references;     /* NULL: none */
  size_t probe_len;
  double probe_error_rate;
  long terminal_window;
  bool exhaustive_search;     /* without the k-mer index (slower) */
  size_t min_length;
  bool tags;
  size_t threads;
} options_t;

/*
   Parse the command line into <self>; the logger is created from the
   common options (-q, -v, -D, -L). Invalid options exit the program
   with an error message and the usage.
*/
void options_parse(options_t *self, int argc, char *argv[],
    Logger **logger);

#endif
