/* (c) 2026 Giorgio Gonnella */

#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "options.h"
#include "argparse-def.h"
#include "amplicon.h"
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
  "            MinKNOW run); the trimmed reads of each input are written\n" \
  "            to <outdir>/<name>.fastq.gz, where <name> is the name of\n" \
  "            the file without its suffix, or that of the directory\n" \
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
  "                          bases, or below 1 a fraction of its length\n" \
  "                          (default: " XSTR(OPTIONS_DEFAULT_MAX_ERRORS) \
  ")\n" \
  "   -w, --window <n>       a primer must lie within <n> bases of its\n" \
  "                          end of the read (default: " \
  XSTR(AMPLICON_DEFAULT_WINDOW) ")\n" \
  "   -m, --min-length <n>   shorter reads after trimming are not\n" \
  "                          written (default: " \
  XSTR(OPTIONS_DEFAULT_MIN_LENGTH) ")\n" \
  "   -j, --threads <n>      threads (default: " \
  XSTR(OPTIONS_DEFAULT_THREADS) "; the output does not\n" \
  "                          depend on it)\n" \
  "   --no-tags              do not append the tags to the read titles\n" \
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
  self->threads = OPTIONS_DEFAULT_THREADS;
  static struct option longopts[] = {
    {"outdir", required_argument, NULL, 'o'},
    {"scheme", required_argument, NULL, 's'},
    {"trim", required_argument, NULL, 'T'},
    {"references", required_argument, NULL, 'r'},
    {"max-errors", required_argument, NULL, 'e'},
    {"window", required_argument, NULL, 'w'},
    {"min-length", required_argument, NULL, 'm'},
    {"threads", required_argument, NULL, 'j'},
    {"no-tags", no_argument, NULL, OPT_NO_TAGS},
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
        if (!parse_double(optarg, &x) || x < 0 ||
            (x >= 1 && x != floor(x)))
          OPT_ERR(c, "a fraction below 1 or a whole number of bases "
                  "expected ('%s' found)", optarg);
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
      case OPT_NO_TAGS:
        self->tags = false;
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
  PARSE_LOG("Errors in a primer: %g\n", self->max_errors);
  PARSE_LOG("Window: %ld\n", self->window);
  PARSE_LOG("Minimum length: %zu\n", self->min_length);
  PARSE_LOG("Threads: %zu\n", self->threads);
}
