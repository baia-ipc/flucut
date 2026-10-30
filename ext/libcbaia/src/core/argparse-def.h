/* (c) 2017-2026 Giorgio Gonnella */

#ifndef ARGPARSE_DEF_H
#define ARGPARSE_DEF_H

#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "logger.h"

/* version of the program, output by the -V option; a program defines it
   before including this header */
#ifndef VERSION
#define VERSION "unknown"
#endif

#define STR(S) #S
#define XSTR(S) STR(S)

#define BOOLSTR(BOOL) ((BOOL) ? "T" : "F")

/*
  Common options and error interface for the CLI scripts

  Usage:
  - a function argsparse() is called from main() for parsing the arguments
    using getopt
  - a macro PURPOSE shall describe the usage of the program
  - PARSE_START is used in the declaration section
  - on errors, use one of the following macros, which format the error,
    output it to stderr and exits the program with code EXIT_FAILURE;
    no newline is necessary in the format string
    - PARSE_ERR(fmt, ...) for generic errors
    - OPT_ERR; ARG_ERR for options and positional args, see below
  - INIT_LOGGER after parsing;
    optionally followed by one or multiple
    PARSE_LOG() for the verbose logging of the option/argument values

  Option reserved:
    -V --> version
    -h --> help
    -q, -v, -D --> verbosity
    -L --> log file (the letter can be changed, see LOGFILE_OPT)

  Put this in the options section of the help, after an empty line:
  COMMON_OPTS_HELP

  Parsing of options (getopt while):
  - COMMON_OPTS_STR append this to the string passed to getopt
  - OPT_ERR(C, fmt, ...) is used for errors concerning the -C option
  - PARSE_INTEGER_OPT_MIN(C, DESTPTR, MINVALUE) checks that the value
          of option C is at least MINVALUE, exits with an error otherwise,
          stores the parse value to *DESTPTR (an integer of any type; exits
          with an error if the value does not fit in it)
  - PARSE_INTEGER_OPT_MIN_MAX(C, DESTPTR, MINVALUE, MAXVALUE) same, but
          with also MAXVALUE check
  - COMMON_OPTS_PARSE_CASES this should be put after the last command-specific
          option parsing case

  Parsing of positional arguments (after getopt while):
  - CHECK_NARGS(EXPECTED) check that there are EXPECTED positional args,
                          exit with an error message otherwise
  - CHECK_MIN_NARGS(EXPECTED) check that there are at least EXPECTED positional
                              args, exit with an error message otherwise
  - ARG_ERR(N, fmt, ...) is used for errors concerning the N-th positional arg
                         (use 0-based numbers)
  - PARSE_INTEGER_ARG_MIN(N, DESTPTR, MINVALUE) checks that the value of the
          N-th positional arg is at least MINVALUE, exits with an error
          otherwise, stores the parse value to *DESTPTR
*/

/* output of the help (PURPOSE) to STREAM; a program can define this macro
   before including this header, e.g. if the help is generated */
#ifndef PRINT_PURPOSE
#define PRINT_PURPOSE(STREAM) fprintf(STREAM, PURPOSE, argv[0])
#endif

/* help and version are output to the standard output, the usage after an
   error to the standard error */
#define HELP_AND_EXIT { \
  PRINT_PURPOSE(stdout);\
  exit(EXIT_SUCCESS); \
}

#define PURPOSE_AND_EXIT { \
  PRINT_PURPOSE(stderr);\
  exit(EXIT_FAILURE); \
}

#define VERSION_AND_EXIT { \
  const char *progname = strrchr(argv[0], '/'); \
  fprintf(stdout, "%s " VERSION "\n", \
          progname == NULL ? argv[0] : progname + 1); \
  exit(EXIT_SUCCESS); \
}

#define LOGGER_OPTS_VAR_DECLARE \
  bool verbose = false,\
       quiet = false,\
       debug = false;\
  FILE* logfile = stderr;

#define PARSE_START \
  LOGGER_OPTS_VAR_DECLARE \

#define PARSE_ERR(FMT, ...) do { \
  lprintf(NULL, LOG_ERROR, FMT "\n\n", ##__VA_ARGS__);\
  PURPOSE_AND_EXIT; \
} while (0)

#define OPT_ERR(C, FMT, ...) \
  PARSE_ERR("Option '-%c': " FMT, C, ##__VA_ARGS__)

#define ARG_ERR(N, FMT, ...) \
  PARSE_ERR("Positional argument %i: " FMT, N, ##__VA_ARGS__)

/* parse a decimal integer (the whole string); false if invalid */
static inline bool argparse_long(const char *s, long *value) {
  char *end;
  errno = 0;
  *value = strtol(s, &end, 10);
  return *s != '\0' && *end == '\0' && errno == 0;
}

/* store VALUE (a long) into *DESTPTR (an integer of any type); false if it
   does not fit, i.e. it is changed by the conversion (also in the sign, for
   an unsigned destination of the size of a long) */
#define ARGPARSE_STORE(DESTPTR, VALUE) \
  ((long)(*(DESTPTR) = (VALUE)) == (VALUE) && \
   (*(DESTPTR) > 0) == ((VALUE) > 0))

#define PARSE_INTEGER_OPT_MIN(C, DESTPTR, MINVALUE) do { \
    long argvalue; \
    if (!argparse_long(optarg, &argvalue) || argvalue < (MINVALUE)) \
      OPT_ERR(C, "argument must be an integer >= %ld (%s found)",\
        (long)(MINVALUE), optarg); \
    if (!ARGPARSE_STORE(DESTPTR, argvalue)) \
      OPT_ERR(C, "argument out of range (%s found)", optarg); \
} while (0)

#define PARSE_INTEGER_OPT_MIN_MAX(C, DESTPTR, MINVALUE, MAXVALUE) do { \
    long argvalue; \
    if (!argparse_long(optarg, &argvalue) || argvalue < (MINVALUE) || \
        argvalue > (MAXVALUE)) \
      OPT_ERR(C, "argument must be an integer >= %ld and <= %ld "\
        "(%s found)", (long)(MINVALUE), (long)(MAXVALUE), optarg); \
    if (!ARGPARSE_STORE(DESTPTR, argvalue)) \
      OPT_ERR(C, "argument out of range (%s found)", optarg); \
} while (0)

#define PARSE_INTEGER_ARG_MIN(N, DESTPTR, MINVALUE) do { \
    long argvalue; \
    if (!argparse_long(argv[optind + (N)], &argvalue) || \
        argvalue < (MINVALUE)) \
      ARG_ERR((N)+1, "argument must be an integer >= %ld (%s found)",\
        (long)(MINVALUE), argv[optind + (N)]); \
    if (!ARGPARSE_STORE(DESTPTR, argvalue)) \
      ARG_ERR((N)+1, "argument out of range (%s found)", \
        argv[optind + (N)]); \
} while (0)

/* option for the log file; a program can define these macros before
   including this header, if it needs the letter for another option */
#ifndef LOGFILE_OPT
#define LOGFILE_OPT 'L'
#define LOGFILE_OPT_STR "L:"
#define LOGFILE_OPT_HELP "-L"
#endif

/*
  Long options: programs using getopt_long() list the common long options
  in their table of options, e.g.:

    static struct option longopts[] = {
      {"threads", required_argument, NULL, 't'},
      COMMON_LONGOPTS
      LONGOPTS_END
    };
    while ((c = getopt_long(argc, argv, COMMON_OPTS_STR "t:", longopts,
                            NULL)) != -1) ...

  (each long option returns the letter of the corresponding short option)
*/
#define LOGGER_LONGOPTS \
  {"quiet", no_argument, NULL, 'q'}, \
  {"verbose", no_argument, NULL, 'v'}, \
  {"debug", no_argument, NULL, 'D'}, \
  {"logfile", required_argument, NULL, LOGFILE_OPT},

#define COMMON_LONGOPTS \
  LOGGER_LONGOPTS \
  {"help", no_argument, NULL, 'h'}, \
  {"version", no_argument, NULL, 'V'},

#define LONGOPTS_END {NULL, 0, NULL, 0}

#define LOGGER_OPTS_HELP \
  "Verbosity options: (only one can be used)\n"\
  "   (default)            success, warning and error messages to stderr\n"\
  "   -q, --quiet          quiet: only error messages\n"\
  "   -v, --verbose        verbose: show additional info messages\n"\
  "   -D, --debug          debug: show some debug information\n"\
  "   " LOGFILE_OPT_HELP ", --logfile F      print messages to file "\
  "(default: stderr)\n"

#define COMMON_OPTS_HELP \
  LOGGER_OPTS_HELP \
  "\n"\
  "Other options:\n"\
  "   -h, --help           show this help message and exit\n"\
  "   -V, --version        output version and exit\n"

#define LOGGER_OPTS_STR "qvD" LOGFILE_OPT_STR

#define COMMON_OPTS_STR LOGGER_OPTS_STR "hV"

#define LOGGER_OPTS_CHKMULTI do { \
  if (verbose || quiet || debug) \
    PARSE_ERR("only one of -v, -q and -D can be used at once"); \
} while (0)

#define LOGGER_OPTS_PARSE_CASES \
  case 'v': \
    LOGGER_OPTS_CHKMULTI; \
    verbose = true; \
    break; \
  case 'q': \
    LOGGER_OPTS_CHKMULTI; \
    quiet = true; \
    break; \
  case 'D': \
    LOGGER_OPTS_CHKMULTI; \
    debug = true; \
    break;\
  case LOGFILE_OPT: \
    /* the option given again: the last one is used */ \
    if (logfile != stderr) \
      fclose(logfile); \
    logfile = fopen(optarg, "w");\
    if (logfile == NULL)\
      OPT_ERR(LOGFILE_OPT, "Opening log file '%s' failed, %s",\
          optarg, strerror(errno));\
    break;

#define COMMON_OPTS_PARSE_CASES \
  LOGGER_OPTS_PARSE_CASES \
  case 'V': \
    VERSION_AND_EXIT; \
    break; \
  case 'h': \
    HELP_AND_EXIT; \
    break; \
  case '?': \
    PARSE_ERR("Unknown option or missing argument: '%s'", \
              argv[optind - 1]); \
    break; \
  default: \
    /* an option of the getopt string without a case */ \
    PARSE_ERR("Option '%s' not handled (internal error)", \
              argv[optind - 1])

#define CHECK_NARGS(EXPECTED) do { \
  if (argc-optind != (EXPECTED)) \
    PARSE_ERR("Wrong number of arguments (%d required, %d found)", \
        (EXPECTED), argc-optind); \
} while (0)

#define CHECK_MIN_NARGS(EXPECTED) do { \
  if (argc-optind < (EXPECTED)) \
    PARSE_ERR("Wrong number of arguments (at least %d required, %d found)", \
        (EXPECTED), argc-optind); \
} while (0)

#define LOGGER_OPTS_APPLY(LOGGER) \
  {\
    LoggingLevel level = quiet ?   LOG_ERROR : \
                         verbose ? LOG_PARAM : \
                         debug ?   LOG_DEBUG : \
                                   LOG_SUCCESS; \
    (LOGGER) = new_logger(logfile, level);\
  }

#define INIT_LOGGER(LOGGER) \
  LOGGER_OPTS_APPLY(LOGGER); \
  Logger *parse_logger = LOGGER; \
  (void)parse_logger;

#define PARSE_LOG(FMT, ...) do { \
  lprintf(parse_logger, LOG_PARAM, FMT, ##__VA_ARGS__); \
} while (0)

#define PARSE_LOG_N(FMT, ...) do { \
  lnprintf(parse_logger, LOG_PARAM, FMT, ##__VA_ARGS__); \
} while (0)

#endif
