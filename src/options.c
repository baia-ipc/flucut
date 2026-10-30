/* (c) 2026 Giorgio Gonnella */

#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include "options.h"
#include "argparse-def.h"
#include "amplicon.h"
#include "auto-settings.h"
#include "live.h"
#include "strutils.h"
#include "worker.h"
#include "reference.h"
#include "scheme.h"

#define OPTIONS_DEFAULT_THREADS 4

/* the long options without a short form */
enum {
  OPT_PROBE_LEN = 256,
  OPT_PROBE_ERRORS,
  OPT_TERMINAL_WINDOW,
  OPT_NO_TAGS,
  OPT_EXHAUSTIVE,
  OPT_ADAPTERS,
  OPT_LIVE,
  OPT_TIMEOUT,
  OPT_MIN_AGE,
  OPT_SCAN_INTERVAL,
  OPT_RESUME,
  OPT_COMPRESS_LEVEL,
};

#define PURPOSE "Remove the primers of the influenza A multisegment " \
  "RT-PCR (Uni12/Uni13)\n" \
  "from the reads, before the assembly (e.g. with IRMA)\n" \
  "\n" \
  "Usage:\n" \
  "   %s [options] <input>...\n" \
  "\n" \
  "Arguments:\n" \
  "   <input>  FASTQ (optionally gzipped) or unaligned BAM file, or a\n" \
  "            directory of such files (e.g. a barcode directory of a\n" \
  "            MinKNOW run); a directory with barcode* subdirectories\n" \
  "            (e.g. fastq_pass) gives one output for each of them; the\n" \
  "            trimmed reads of each are written to\n" \
  "            <outdir>/<name>.fastq.gz, where <name> is the name of the\n" \
  "            file without its suffix, or that of the directory\n" \
  "\n" \
  "Options:\n" \
  "   -o, --outdir <dir>     output directory (default: " \
  OPTIONS_DEFAULT_OUTDIR ")\n" \
  "   -s, --scheme <name>    primer scheme (default: " SCHEME_DEFAULT \
  "; see\n" \
  "                          below)\n" \
  "   -T, --trim <policy>    primer: remove every primer-derived base\n" \
  "                          (default); terminus: keep the conserved\n" \
  "                          termini (requires --references)\n" \
  "   -r, --references <f>   FASTA of segment references: the reads are\n" \
  "                          classified as terminal or internally primed\n" \
  "                          (with --trim primer: reported only)\n" \
  "   -e, --max-errors <x>   errors allowed in a primer: a number of\n" \
  "                          bases, below 1 a fraction of its length, or\n" \
  "                          auto: " XSTR(AUTO_MAX_ERRORS_LOW) ", or " \
  XSTR(AUTO_MAX_ERRORS_HIGH) " if the primers of the\n" \
  "                          first reads have at least " \
  XSTR(AUTO_MAX_ERRORS_HIGH_RATE) " errors per\n" \
  "                          base (default: auto)\n" \
  "   -w, --window <n>       a primer must lie within <n> bases of its\n" \
  "                          end of the read (default: " \
  XSTR(AMPLICON_DEFAULT_WINDOW) ")\n" \
  "   -m, --min-length <n>   shorter reads after trimming are not\n" \
  "                          written (default: " \
  XSTR(OPTIONS_DEFAULT_MIN_LENGTH) ")\n" \
  "   -j, --threads <n>      threads (default: " \
  XSTR(OPTIONS_DEFAULT_THREADS) "; the output does not\n" \
  "                          depend on it)\n" \
  "   --compress-level <n>   gzip compression level of the outputs, " \
  XSTR(WORKER_MIN_COMPRESS_LEVEL) " to\n" \
  "                          " XSTR(WORKER_MAX_COMPRESS_LEVEL) \
  " (default: " XSTR(WORKER_DEFAULT_COMPRESS_LEVEL) "; higher: smaller, " \
  "slower)\n" \
  "   --no-tags              do not append the tags to the read titles\n" \
  "   --adapters <kits>      the Nanopore adapters searched at the read\n" \
  "                          ends without a primer: a comma-separated\n" \
  "                          list of " ADAPTERS_KITS ", or none\n" \
  "                          (default: " ADAPTERS_DEFAULT ")\n" \
  "\n" \
  "Live mode (a MinKNOW run in progress):\n" \
  "   --live                 process the files of the input directories as\n" \
  "                          they are written, appending to the outputs and\n" \
  "                          updating the tables and reports after each\n" \
  "                          file; stop at the end of the run (a file\n" \
  "                          " LIVE_END_FILES " in an input\n" \
  "                          directory or its parent), after --timeout, or\n" \
  "                          with Ctrl-C\n" \
  "   --timeout <t>          also stop after <t> without new files (e.g.\n" \
  "                          30s, 60m, 2h; a number: minutes; default:\n" \
  "                          none)\n" \
  "   --min-age <t>          process a file only <t> after its last\n" \
  "                          change (e.g. 30s, 2m; a number: seconds;\n" \
  "                          default: " XSTR(OPTIONS_DEFAULT_MIN_AGE) "s)\n" \
  "   --scan-interval <s>    seconds between the scans of the input\n" \
  "                          directories (default: " \
  XSTR(OPTIONS_DEFAULT_SCAN_INTERVAL) ")\n" \
  "   --resume               skip the files already processed by an\n" \
  "                          interrupted live run (same output directory)\n" \
  "\n" \
  "Classification against the references (with --references):\n" \
  "   --probe-len <n>        bases after the primer located in the\n" \
  "                          references (default: " \
  XSTR(REFERENCE_DEFAULT_PROBE_LEN) ", minimum: " \
  XSTR(REFERENCE_MIN_PROBE_LEN) ")\n" \
  "   --probe-errors <x>     fraction of the probe allowed to be errors\n" \
  "                          (default: " \
  XSTR(REFERENCE_DEFAULT_MAX_ERROR_RATE) ")\n" \
  "   --terminal-window <n>  a probe within <n> bases of a segment end\n" \
  "                          is terminal (default: " \
  XSTR(REFERENCE_DEFAULT_TERMINAL_WINDOW) ")\n" \
  "   --exhaustive-search    align each probe to the whole of every\n" \
  "                          reference, without the k-mer index (slower;\n" \
  "                          for checking the index)\n" \
  "\n" \
  COMMON_OPTS_HELP

/* the help, with the list of the built-in schemes */
static void print_purpose(FILE *stream, const char *progname) {
  fprintf(stream, PURPOSE, progname);
  fputs("\nPrimer schemes:", stream);
  for (const char *const *n = scheme_builtin_names(); *n != NULL; n++)
    fprintf(stream, " %s", *n);
  fputs("\n", stream);
}

#undef PRINT_PURPOSE
#define PRINT_PURPOSE(STREAM) print_purpose(STREAM, argv[0])

/* parse a decimal number (the whole string); false if invalid */
static bool parse_double(const char *s, double *value) {
  char *end;
  errno = 0;
  *value = strtod(s, &end);
  return *s != '\0' && *end == '\0' && errno == 0 && isfinite(*value);
}

/* whether <name> is a built-in scheme (ignoring the case) */
static bool is_scheme(const char *name) {
  for (const char *const *n = scheme_builtin_names(); *n != NULL; n++)
    if (strcasecmp(*n, name) == 0)
      return true;
  return false;
}

void options_parse(options_t *self, int argc, char *argv[],
    Logger **logger) {
  PARSE_START;
  self->outdir = OPTIONS_DEFAULT_OUTDIR;
  self->scheme = SCHEME_DEFAULT;
  self->policy = AMPLICON_DEFAULT_POLICY;
  self->max_errors = OPTIONS_DEFAULT_MAX_ERRORS;
  self->window = AMPLICON_DEFAULT_WINDOW;
  self->references = NULL;
  self->probe_len = REFERENCE_DEFAULT_PROBE_LEN;
  self->probe_error_rate = REFERENCE_DEFAULT_MAX_ERROR_RATE;
  self->terminal_window = REFERENCE_DEFAULT_TERMINAL_WINDOW;
  self->exhaustive_search = false;
  self->min_length = OPTIONS_DEFAULT_MIN_LENGTH;
  self->tags = true;
  self->adapters = ADAPTERS_DEFAULT;
  self->live = false;
  self->live_timeout = 0;
  self->live_min_age = OPTIONS_DEFAULT_MIN_AGE;
  self->live_scan_interval = OPTIONS_DEFAULT_SCAN_INTERVAL;
  self->resume = false;
  bool timeout_given = false, min_age_given = false,
       scan_interval_given = false;
  self->threads = OPTIONS_DEFAULT_THREADS;
  self->compress_level = WORKER_DEFAULT_COMPRESS_LEVEL;
  static struct option longopts[] = {
    {"outdir", required_argument, NULL, 'o'},
    {"scheme", required_argument, NULL, 's'},
    {"trim", required_argument, NULL, 'T'},
    {"references", required_argument, NULL, 'r'},
    {"max-errors", required_argument, NULL, 'e'},
    {"window", required_argument, NULL, 'w'},
    {"min-length", required_argument, NULL, 'm'},
    {"threads", required_argument, NULL, 'j'},
    {"compress-level", required_argument, NULL, OPT_COMPRESS_LEVEL},
    {"no-tags", no_argument, NULL, OPT_NO_TAGS},
    {"adapters", required_argument, NULL, OPT_ADAPTERS},
    {"live", no_argument, NULL, OPT_LIVE},
    {"timeout", required_argument, NULL, OPT_TIMEOUT},
    {"min-age", required_argument, NULL, OPT_MIN_AGE},
    {"scan-interval", required_argument, NULL, OPT_SCAN_INTERVAL},
    {"resume", no_argument, NULL, OPT_RESUME},
    {"probe-len", required_argument, NULL, OPT_PROBE_LEN},
    {"probe-errors", required_argument, NULL, OPT_PROBE_ERRORS},
    {"terminal-window", required_argument, NULL, OPT_TERMINAL_WINDOW},
    {"exhaustive-search", no_argument, NULL, OPT_EXHAUSTIVE},
    COMMON_LONGOPTS
    LONGOPTS_END
  };
  int c;
  double x;
  opterr = 0;
  while ((c = getopt_long(argc, argv, COMMON_OPTS_STR "o:s:T:r:e:w:m:j:",
                          longopts, NULL)) != -1) {
    switch (c) {
      case 'o':
        if (optarg[0] == '\0')
          OPT_ERR(c, "the output directory is empty");
        self->outdir = optarg;
        break;
      case 's':
        if (!is_scheme(optarg))
          OPT_ERR(c, "unknown scheme '%s' (see the list below)", optarg);
        self->scheme = optarg;
        break;
      case 'T':
        if (strcmp(optarg, "primer") == 0)
          self->policy = TRIM_PRIMER;
        else if (strcmp(optarg, "terminus") == 0)
          self->policy = TRIM_TERMINUS;
        else
          OPT_ERR(c, "'primer' or 'terminus' expected ('%s' found)",
                  optarg);
        break;
      case 'r':
        self->references = optarg;
        break;
      case 'e':
        if (strcmp(optarg, "auto") == 0) {
          self->max_errors = OPTIONS_MAX_ERRORS_AUTO;
          break;
        }
        if (!parse_double(optarg, &x) || x < 0 ||
            (x >= 1 && x != floor(x)))
          OPT_ERR(c, "a fraction below 1, a whole number of bases or "
                  "auto expected ('%s' found)", optarg);
        self->max_errors = x;
        break;
      case 'w':
        PARSE_INTEGER_OPT_MIN(c, &self->window, 0);
        break;
      case 'm':
        PARSE_INTEGER_OPT_MIN(c, &self->min_length, 1);
        break;
      case 'j':
        PARSE_INTEGER_OPT_MIN_MAX(c, &self->threads, 1,
                                  OPTIONS_MAX_THREADS);
        break;
      case OPT_COMPRESS_LEVEL:
        {
          long v;
          if (!argparse_long(optarg, &v) || v < WORKER_MIN_COMPRESS_LEVEL ||
              v > WORKER_MAX_COMPRESS_LEVEL)
            PARSE_ERR("Option '--compress-level': an integer >= "
                      XSTR(WORKER_MIN_COMPRESS_LEVEL) " and <= "
                      XSTR(WORKER_MAX_COMPRESS_LEVEL) " expected ('%s' "
                      "found)", optarg);
          self->compress_level = (int)v;
        }
        break;
      case OPT_NO_TAGS:
        self->tags = false;
        break;
      case OPT_LIVE:
        self->live = true;
        break;
      case OPT_TIMEOUT:
        if (!parse_duration(optarg, 60, &self->live_timeout) ||
            self->live_timeout == 0)
          PARSE_ERR("Option '--timeout': a positive duration expected, "
                    "e.g. 30s, 60m, 2h ('%s' found)", optarg);
        timeout_given = true;
        break;
      case OPT_MIN_AGE:
        if (!parse_duration(optarg, 1, &self->live_min_age))
          PARSE_ERR("Option '--min-age': a duration expected, e.g. 30s, "
                    "2m ('%s' found)", optarg);
        min_age_given = true;
        break;
      case OPT_SCAN_INTERVAL:
        if (!parse_double(optarg, &x) || x <= 0 ||
            x > OPTIONS_MAX_SCAN_INTERVAL)
          PARSE_ERR("Option '--scan-interval': a number of seconds > 0 and "
                    "<= " XSTR(OPTIONS_MAX_SCAN_INTERVAL) " expected ('%s' "
                    "found)", optarg);
        self->live_scan_interval = x;
        scan_interval_given = true;
        break;
      case OPT_RESUME:
        self->resume = true;
        break;
      case OPT_ADAPTERS:
        if (!adapters_valid(optarg))
          PARSE_ERR("Option '--adapters': a comma-separated list of "
                    ADAPTERS_KITS ", or none, expected ('%s' found)",
                    optarg);
        self->adapters = optarg;
        break;
      case OPT_EXHAUSTIVE:
        self->exhaustive_search = true;
        break;
      case OPT_PROBE_LEN:
        {
          long v;
          if (!argparse_long(optarg, &v) || v < REFERENCE_MIN_PROBE_LEN)
            PARSE_ERR("Option '--probe-len': an integer >= %d expected "
                      "('%s' found)", REFERENCE_MIN_PROBE_LEN, optarg);
          self->probe_len = (size_t)v;
        }
        break;
      case OPT_PROBE_ERRORS:
        if (!parse_double(optarg, &x) || x < 0 || x >= 1)
          PARSE_ERR("Option '--probe-errors': a fraction >= 0 and < 1 "
                    "expected ('%s' found)", optarg);
        self->probe_error_rate = x;
        break;
      case OPT_TERMINAL_WINDOW:
        {
          long v;
          if (!argparse_long(optarg, &v) || v < 0)
            PARSE_ERR("Option '--terminal-window': an integer >= 0 "
                      "expected ('%s' found)", optarg);
          self->terminal_window = v;
        }
        break;
      COMMON_OPTS_PARSE_CASES;
    }
  }
  CHECK_MIN_NARGS(1);
  if (self->policy == TRIM_TERMINUS && self->references == NULL)
    PARSE_ERR("--trim terminus requires --references: without them the "
              "internally primed reads cannot be recognised, and would "
              "keep the primer bases (the artefact)");
  /* the options of the live mode require it */
  if (timeout_given && !self->live)
    PARSE_ERR("Option '--timeout': the live mode (--live) is required");
  if (min_age_given && !self->live)
    PARSE_ERR("Option '--min-age': the live mode (--live) is required");
  if (scan_interval_given && !self->live)
    PARSE_ERR("Option '--scan-interval': the live mode (--live) is "
              "required");
  if (self->resume && !self->live)
    PARSE_ERR("Option '--resume': the live mode (--live) is required");
  if (self->live)
    for (int i = optind; i < argc; i++) {
      struct stat st;
      if (stat(argv[i], &st) != 0 || !S_ISDIR(st.st_mode))
        PARSE_ERR("Option '--live': the inputs must be directories ('%s' "
                  "found)", argv[i]);
    }
  self->inputs = argv + optind;
  self->n_inputs = (size_t)(argc - optind);
  self->show_banner = !quiet && logfile == stderr &&
                      isatty(fileno(stderr));
  INIT_LOGGER(*logger);
  PARSE_LOG("Inputs: %zu\n", self->n_inputs);
  PARSE_LOG("Output directory: %s\n", self->outdir);
  PARSE_LOG("Scheme: %s\n", self->scheme);
  PARSE_LOG("Trimming: %s\n",
            self->policy == TRIM_PRIMER ? "primer" : "terminus");
  PARSE_LOG("References: %s\n",
            self->references == NULL ? "none" : self->references);
  if (self->max_errors == OPTIONS_MAX_ERRORS_AUTO)
    PARSE_LOG("Errors in a primer: auto\n");
  else
    PARSE_LOG("Errors in a primer: %g\n", self->max_errors);
  PARSE_LOG("Window: %ld\n", self->window);
  PARSE_LOG("Minimum length: %zu\n", self->min_length);
  PARSE_LOG("Adapters: %s\n", self->adapters);
  if (self->live) {
    PARSE_LOG("Live mode: timeout %ld s (0: none), minimum age %ld s, "
              "scan interval %g s%s\n", self->live_timeout,
              self->live_min_age, self->live_scan_interval,
              self->resume ? ", resume" : "");
  }
  PARSE_LOG("Threads: %zu\n", self->threads);
  PARSE_LOG("Compression level: %d\n", self->compress_level);
}
