/* (c) 2026 Giorgio Gonnella */

/*
   Example program of the libcbaia manual: the messages of the logger at
   each level, and the progress display.
*/

#define VERSION "1.0"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "argparse-def.h"
#include "logger.h"
#include "progress.h"

#define DEFAULT_N 2000000

// %-PURPOSE-%
#define PURPOSE "Show the messages of the logger and a progress display\n"\
  "\n"\
  "Usage:\n"\
  "   %s [options]\n"\
  "\n"\
  "Options:\n"\
  "   -n, --items <n>  number of items to process (default: "\
  XSTR(DEFAULT_N) ")\n"\
  "   -p, --progress   show a progress display\n"\
  "   -e, --error      end with an error\n"\
  "\n"\
  COMMON_OPTS_HELP

void parseargs(int argc, char *argv[],
    long *n_items,
    bool *show_progress,
    bool *error,
    Logger **logger)
{
  PARSE_START;
  int c = 0;
  *n_items = DEFAULT_N;
  *show_progress = false;
  *error = false;
  static struct option longopts[] = {
    {"items", required_argument, NULL, 'n'},
    {"progress", no_argument, NULL, 'p'},
    {"error", no_argument, NULL, 'e'},
    COMMON_LONGOPTS
    LONGOPTS_END
  };
  opterr = 0;
  while ((c = getopt_long(argc, argv, COMMON_OPTS_STR "n:pe", longopts,
                          NULL)) != -1) {
    switch (c) {
      case 'n':
        PARSE_INTEGER_OPT_MIN('n', n_items, 1);
        break;
      case 'p':
        *show_progress = true;
        break;
      case 'e':
        *error = true;
        break;
      COMMON_OPTS_PARSE_CASES;
    }
  }
  CHECK_NARGS(0);
  INIT_LOGGER(*logger);
  PARSE_LOG("Items: %ld\n", *n_items);
  PARSE_LOG("Progress display: %s\n", BOOLSTR(*show_progress));
}

int main(int argc, char *argv[]) {
  long n_items;
  bool show_progress, error;
  Logger *logger;
  parseargs(argc, argv, &n_items, &show_progress, &error, &logger);
  lprintf(logger, LOG_PHASE, "Processing the items...\n");
  progress_t *progress = show_progress ?
    progress_new(stderr, "Processing items", (size_t)n_items) : NULL;
  size_t n_odd = 0;
  for (long i = 0; i < n_items; i++) {
    if (i % 7 == 1)
      n_odd++;
    if (progress != NULL && i % 1000 == 999) {
      usleep(500);   /* (the work of the example is too fast to see) */
      progress_update(progress, 1000);
    }
  }
  if (progress != NULL) {
    progress_update(progress, (size_t)(n_items % 1000));
    progress_finish(progress);
    progress_delete(progress);
  }
  lprintf(logger, LOG_INFO, "Items with remainder 1 (mod 7): %zu\n",
          n_odd);
  lprintf(logger, LOG_DEBUG, "Last item: %ld\n", n_items - 1);
  showmem(logger);
  lprintf(logger, LOG_WARN, "%zu items were skipped\n", (size_t)3);
  LCHECK(!error, logger, "Item %ld: invalid value (-e was given)\n",
         n_items / 2);
  lprintf(logger, LOG_SUCCESS, "%ld items processed\n", n_items);
  delete_logger(logger);
  return EXIT_SUCCESS;
}
