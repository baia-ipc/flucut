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
/* -e auto (see auto-settings.h): the default */
#define OPTIONS_MAX_ERRORS_AUTO -1.0
#define OPTIONS_DEFAULT_MAX_ERRORS OPTIONS_MAX_ERRORS_AUTO
#define OPTIONS_DEFAULT_MIN_LENGTH 1
#define OPTIONS_MAX_THREADS 256
/* --platform: the settings of each platform (see the help): the window
   of the Illumina reads (in which the primer is at the read start, the
   adapters removed by the instrument) */
#define OPTIONS_ILLUMINA_WINDOW 10
/* the live mode (see live.h): defaults of --min-age (seconds) and
   --scan-interval (seconds), maximum of --scan-interval */
#define OPTIONS_DEFAULT_MIN_AGE 60
#define OPTIONS_DEFAULT_SCAN_INTERVAL 10
#define OPTIONS_MAX_SCAN_INTERVAL 3600

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
  /* --platform: illumina (false: nanopore) */
  bool illumina;
  /* the user's references file (-r <file>); NULL and references_off false:
     the built-in set (the default, see builtin-references.h) */
  const char *references;
  bool references_off;        /* -r off: no references (primer-only) */
  size_t probe_len;
  double probe_error_rate;
  long terminal_window;
  bool exhaustive_search;     /* without the k-mer index (slower) */
  size_t min_length;
  bool tags;
  const char *adapters;       /* see adapters.h */
  size_t threads;
  int compress_level;         /* of the outputs (gzip) */
  /* the live mode (see live.h): --live, --timeout (seconds, 0: none),
     --min-age (seconds), --scan-interval (seconds), --resume */
  bool live;
  long live_timeout;
  long live_min_age;
  double live_scan_interval;
  bool resume;
  /* show the banner: not quiet, messages to the standard error, which
     is a terminal */
  bool show_banner;
} options_t;

/*
   Parse the command line into <self>; the logger is created from the
   common options (-q, -v, -D, -L). Invalid options exit the program
   with an error message and the usage.
*/
void options_parse(options_t *self, int argc, char *argv[],
    Logger **logger);

#endif
