/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <dirent.h>
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
/* first: the help (with the test options of the table of the test
   groups, see print_purpose()) and the version (VERSION, options.h), used
   by argparse-def.h */
#define PRINT_PURPOSE(STREAM) print_purpose(STREAM, argv[0])
#include "options.h"
#include "argparse-def.h"
#include "adapters.h"
#include "adaptor-search.h"
#include "amplicon.h"
#include "amplicon-test.h"
#include "anchor-site.h"
#include "auto-settings.h"
#include "platform-check.h"
#include "qc-flags.h"
#include "qc-summary.h"
#include "run-resources.h"
#include "read-metadata.h"
#include "read-stats.h"
#include "chimera.h"
#include "cleavage-counts-test.h"
#include "codons.h"
#include "counter.h"
#include "fasta.h"
#include "lengths.h"
#include "file.h"
#include "input-groups.h"
#include "internal-matches.h"
#include "flank-deletions.h"
#include "junction.h"
#include "logger.h"
#include "markers-test.h"
#include "placement.h"
#include "read-agreement.h"
#include "virus.h"
#include "read-cover.h"
#include "read-outputs.h"
#include "primer-search.h"
#include "primer-search-test.h"
#include "primers.h"
#include "reference.h"
#include "reference-test.h"
#include "report.h"
#include "scheme.h"
#include "scheme-choice.h"
#include "segments.h"
#include "subtype-reads.h"
#include "terminus-matches.h"
#include "trimmer.h"
#include "user-config.h"

/* the directory for temporary files (-T; NULL: a new one, see
   make_tmpdir()) */
static const char *tmpdir = NULL;

static void run_input_groups_test(void) { input_groups_test(tmpdir); }
static void run_read_outputs_test(void) { read_outputs_test(tmpdir); }

/* a new directory for the temporary files, in $TMPDIR or /tmp (so that
   several runs do not use the same files); <path>: its path */
static void make_tmpdir(char *path, size_t size, Logger *logger) {
  const char *base = getenv("TMPDIR");
  if (base == NULL || base[0] == '\0')
    base = "/tmp";
  bool ok = (size_t)snprintf(path, size, "%s/fc-tests-XXXXXX", base) < size
            && mkdtemp(path) != NULL;
  LCHECK(ok, logger, "Cannot create a temporary directory in %s: %s\n",
         base, strerror(errno));
}

/* remove the directory <path> and the files in it (the tests write no
   subdirectories) */
static void remove_tmpdir(const char *path, Logger *logger) {
  DIR *dir = opendir(path);
  bool ok = dir != NULL;
  for (struct dirent *e = ok ? readdir(dir) : NULL; e != NULL;
       e = readdir(dir)) {
    if (strcmp(e->d_name, ".") == 0 || strcmp(e->d_name, "..") == 0)
      continue;
    char *file = path_join(path, e->d_name);
    ok = remove(file) == 0 && ok;
    free(file);
  }
  if (dir != NULL)
    closedir(dir);
  if (!ok || rmdir(path) != 0)
    lprintf(logger, LOG_WARN, "Cannot remove the temporary directory %s: "
            "%s\n", path, strerror(errno));
}

typedef struct {
  char option;
  const char *description;
  void (*run)(void);
  bool selected;
} test_group_t;

static test_group_t test_groups[] = {
  {'p', "the built-in primers", primers_test, false},
  {'s', "the search of a primer in a read", primer_search_test, false},
  {'i', "the matches of the primers inside a read",
   internal_matches_test, false},
  {'c', "the primer schemes", scheme_test, false},
  {'a', "the analysis and trimming of a read", amplicon_test, false},
  {'r', "the location of the reads in the references", reference_test,
   false},
  {'P', "the placement of a read on the references", placement_test,
   false},
  {'J', "the junction of a read with an internal deletion",
   junction_test, false},
  {'F', "the short deletions next to the junction",
   flank_deletions_test, false},
  {'O', "the category of a read, from where it lies on its reference",
   read_cover_test, false},
  {'e', "the matches of the conserved termini in a read",
   terminus_matches_test, false},
  {'C', "the splitting of the reads with several amplicons",
   chimera_test, false},
  {'t', "the trimming of a batch of reads", trimmer_test, false},
  {'l', "the lengths of the complete amplicons", lengths_test, false},
  {'n', "the counts of the reads", counter_test, false},
  {'R', "the content of the reports", report_test, false},
  {'g', "the groups of input files", run_input_groups_test, false},
  {'o', "the output files of a group", run_read_outputs_test, false},
  {'S', "the recognition of the segments", segments_test, false},
  {'B', "the type of the virus, A or B", virus_test, false},
  {'k', "the genetic code", codons_test, false},
  {'E', "the anchor of a site inside a read", anchor_site_test, false},
  {'H', "the cleavage site of the haemagglutinin in the reads",
   cleavage_counts_test, false},
  {'A', "the search of the adapters", adapters_test, false},
  {'u', "the automatic settings", auto_settings_test, false},
  {'W', "the platform read from the reads", platform_check_test, false},
  {'x', "the choice of the primer scheme from the reads",
   scheme_choice_test, false},
  {'G', "the agreement of the probes of a read", read_agreement_test,
   false},
  {'b', "the reads of a second HA or NA subtype", subtype_reads_test,
   false},
  {'Q', "the summary of a sample by segment", qc_summary_test, false},
  {'f', "what deserves attention in a sample", qc_flags_test, false},
  {'d', "the search of an adaptor", adaptor_search_test,
   false},
  {'U', "the resources used by a run", run_resources_test, false},
  {'Y', "the configuration file", user_config_test, false},
  {'m', "the markers read from the reads", markers_test, false},
};
#define N_TEST_GROUPS (sizeof(test_groups) / sizeof(test_groups[0]))

static void print_purpose(FILE *stream, const char *progname) {
  fprintf(stream, "Unit tests of flucut\n"
          "\n"
          "Usage:\n"
          "   %s [options]\n"
          "\n"
          "Test options (default: run all tests):\n", progname);
  for (size_t i = 0; i < N_TEST_GROUPS; i++)
    fprintf(stream, " -%c : test %s\n", test_groups[i].option,
            test_groups[i].description);
  fputs("\n"
        "Other options:\n"
        " -T <dir> : directory for temporary files (default: a new directory\n"
        "            in $TMPDIR, or /tmp, removed at the end)\n"
        "\n"
        COMMON_OPTS_HELP, stream);
}

void parseargs(int argc, char *argv[], Logger **logger) {
  PARSE_START;
  int c = 0;
  bool any = false;
  opterr = 0;
  /* the other options, then one letter for each test group */
  char optstring[sizeof(COMMON_OPTS_STR "T:") + N_TEST_GROUPS] =
    COMMON_OPTS_STR "T:";
  size_t len = strlen(optstring);
  for (size_t i = 0; i < N_TEST_GROUPS; i++) {
    /* The letter of a test group must be its own: not one of the
       options every program has (-q, -v, -D, -L, -h, -V, -T), and not
       one already taken by another group, or it would select the wrong
       group or none. Checked here and not with assert(), so that it is
       checked in every build and a table which cannot work is a plain
       error rather than a silent surprise. */
    if (strchr(optstring, test_groups[i].option) != NULL) {
      fprintf(stderr, "%s: the letter '%c' of the test group \"%s\" is "
              "already an option\n", argv[0], test_groups[i].option,
              test_groups[i].description);
      exit(EXIT_FAILURE);
    }
    assert(len + 1 < sizeof(optstring));
    optstring[len++] = test_groups[i].option;
    optstring[len] = '\0';
  }
  while ((c = getopt(argc, argv, optstring)) != -1) {
    bool found = false;
    for (size_t i = 0; i < N_TEST_GROUPS; i++)
      if (c == test_groups[i].option) {
        test_groups[i].selected = true;
        any = found = true;
      }
    if (found)
      continue;
    switch (c) {
      case 'T':
        tmpdir = optarg;
        break;
      COMMON_OPTS_PARSE_CASES;
    }
  }
  CHECK_NARGS(0);
  if (!any)
    for (size_t i = 0; i < N_TEST_GROUPS; i++)
      test_groups[i].selected = true;
  INIT_LOGGER(*logger);
}

int main(int argc, char *argv[]) {
  Logger *logger;
  parseargs(argc, argv, &logger);
  char own_tmpdir[4096];
  bool own = tmpdir == NULL;
  if (own) {
    make_tmpdir(own_tmpdir, sizeof(own_tmpdir), logger);
    tmpdir = own_tmpdir;
  }
  for (size_t i = 0; i < N_TEST_GROUPS; i++)
    if (test_groups[i].selected) {
      lprintf(logger, LOG_PHASE, "Testing %s...\n",
              test_groups[i].description);
      test_groups[i].run();
    }
  if (own)
    remove_tmpdir(own_tmpdir, logger);
  lprintf(logger, LOG_SUCCESS, "All tests passed\n");
  delete_logger(logger);
  return EXIT_SUCCESS;
}
