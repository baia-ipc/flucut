/* (c) 2026 Giorgio Gonnella */

/* version of libcbaia */
#define VERSION "0.1.0"

/* the help, with the test options of the table of the test groups */
#define PRINT_PURPOSE(STREAM) print_purpose(STREAM, argv[0])

#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "argparse-def.h"
#include "bam-reader.h"
#include "bgzf-compressor.h"
#include "bitwise.h"
#include "bloom-filter.h"
#include "cell-caller.h"
#include "complexity.h"
#include "fastq-reader.h"
#include "fastq-writer.h"
#include "file.h"
#include "file-merger.h"
#include "file-sorter.h"
#include "gzip-compressor.h"
#include "heap.h"
#include "html-page.h"
#include "html-template.h"
#include "html-writer.h"
#include "hyperloglog.h"
#include "logger.h"
#include "md5.h"
#include "outfiles.h"
#include "progress.h"
#include "reads-reader.h"
#include "edit-distance.h"
#include "minknow-metadata.h"
#include "numutils.h"
#include "sam-header.h"
#include "sam-tags.h"
#include "sequtil.h"
#include "set.h"
#include "strbuf.h"
#include "strutils.h"
#include "test.h"
#include "ubam-writer.h"
#include "whitelist.h"

#define DEFAULT_TMPDIR "/tmp"

/* the directory for temporary files (option -T), passed to the tests
   which need it */
static const char *opt_tmpdir = DEFAULT_TMPDIR;

static void strutils_test(void) {
  char *s = strcat_alloc("abc", "de");
  test_eq_str(s, "abcde");
  free(s);
  s = xstrdup("abc");
  test_eq_str(s, "abc");
  free(s);
  strbuf_t buf;
  strbuf_init(&buf);
  strbuf_append_str(&buf, "x=");
  strbuf_append_long(&buf, -1234);
  strbuf_append_char(&buf, ';');
  strbuf_append_long(&buf, 0);
  test_assert(buf.len == 9 && memcmp(buf.s, "x=-1234;0", 9) == 0);
  strbuf_truncate(&buf, 7);
  strbuf_append_char(&buf, ',');
  test_assert(buf.len == 8 && memcmp(buf.s, "x=-1234,", 8) == 0);
  strbuf_clear(&buf);
  for (size_t i = 0; i < 10000; i++)
    strbuf_append(&buf, "abc", 3);
  test_eq_size(buf.len, 30000);
  strbuf_free(&buf);
  test_assert(str_has_suffix("reads.fastq.gz", ".gz"));
  test_assert(!str_has_suffix(".gz", ".gz"));
  test_assert(!str_has_suffix("reads.fastq", ".gz"));
  char dur[32];
  format_duration(3725, dur, sizeof(dur));
  test_eq_str(dur, "01:02:05");
  format_duration(90061, dur, sizeof(dur));
  test_eq_str(dur, "1 d 01:01:01");
  format_duration(-3, dur, sizeof(dur));
  test_eq_str(dur, "00:00:00");
  long sec;
  test_assert(parse_duration("30s", 60, &sec) && sec == 30);
  test_assert(parse_duration("2m", 1, &sec) && sec == 120);
  test_assert(parse_duration("1.5h", 1, &sec) && sec == 5400);
  test_assert(parse_duration("5", 60, &sec) && sec == 300);
  test_assert(parse_duration("0", 1, &sec) && sec == 0);
  test_assert(!parse_duration("", 1, &sec));
  test_assert(!parse_duration("-1s", 1, &sec));
  test_assert(!parse_duration("5x", 1, &sec));
  test_assert(!parse_duration("5ms", 1, &sec));
  test_assert(!parse_duration("s", 1, &sec));
  test_assert(!parse_duration("inf", 1, &sec));
  test_assert(!parse_duration("nan", 1, &sec));
  test_assert(!parse_duration("1e30h", 1, &sec));
  test_assert(is_space(' ') && is_space('\t') && is_space('\r') &&
              is_space('\x1f'));
  test_assert(!is_space('a') && !is_space('\0') && !is_space('\x1b'));
  char *out = NULL;
  size_t outsize = 0;
  FILE *f = open_memstream(&out, &outsize);
  test_assert(f != NULL);
  fputs_json_string("a\"b\\c\n\x01", f);
  fputc(' ', f);
  fputs_html_escaped("<a href=\"x\">&</a>", f);
  fclose(f);
  test_eq_str(out, "\"a\\\"b\\\\c\\n\\u0001\" "
              "&lt;a href=&quot;x&quot;&gt;&amp;&lt;/a&gt;");
  free(out);
}

static bool only_txt(const char *name, const char *path, void *data) {
  (void)path;
  (*(size_t *)data)++;
  return str_has_suffix(name, ".txt");
}

static void file_utils_test(const char *tmpdir) {
  char *joined = path_join("out", "S01.fastq.gz");
  test_eq_str(joined, "out/S01.fastq.gz");
  free(joined);
  joined = path_join("out/", "S01");
  test_eq_str(joined, "out/S01");
  free(joined);
  joined = path_join("", "S01");
  test_eq_str(joined, "S01");
  free(joined);
  char *base = path_basename("run/fastq_pass/");
  test_eq_str(base, "fastq_pass");
  free(base);
  base = path_basename("reads.fastq");
  test_eq_str(base, "reads.fastq");
  free(base);
  char *dir = strcat_alloc(tmpdir, "/cbaia-tests-file-utils");
  char *sub = strcat_alloc(dir, "/a/b");
  xmakedirs(sub, NULL);
  xmakedirs(sub, NULL);   /* existing: no error */
  test_assert(is_directory(sub));
  test_assert(!is_directory("/nonexistent/x"));
  const char *names[] = {"z.txt", "y.dat", "x.txt", ".hidden.txt"};
  for (size_t i = 0; i < 4; i++) {
    char *dir_slash = strcat_alloc(dir, "/");
    char *path = strcat_alloc(dir_slash, names[i]);
    FILE *f = fopen(path, "w");
    test_assert(f != NULL);
    fclose(f);
    free(dir_slash);
    free(path);
  }
  /* read_file: an empty file (z.txt), then with a content */
  char *zpath = path_join(dir, "z.txt");
  FILE *zf = fopen(zpath, "r+");
  test_assert(zf != NULL);
  char *content;
  size_t content_size;
  read_file(zf, &content, &content_size, NULL);
  test_eq_size(content_size, 0);
  test_eq_str(content, "");
  free(content);
  test_assert(fputs("abc\n", zf) >= 0);
  test_eq_size(file_size(zf), 4);
  read_file(zf, &content, &content_size, NULL);
  test_eq_size(content_size, 4);
  test_eq_str(content, "abc\n");
  free(content);
  /* read_lines_in_buffer with the minimum buffer size: one empty line
     at a time */
  test_assert(freopen(NULL, "w+", zf) != NULL);
  test_assert(fputs("\n\n", zf) >= 0);
  rewind(zf);
  char buffer[READ_LINES_MIN_BUFSIZE], **lines = NULL;
  size_t linesbufsize = 0, nlines = 0, nread;
  while ((nread = read_lines_in_buffer(zf, buffer, sizeof buffer, &lines,
          &linesbufsize, NULL)) > 0) {
    test_eq_size(nread, 1);
    test_eq_str(lines[0], "\n");
    nlines += nread;
  }
  test_eq_size(nlines, 2);
  free(lines);
  fclose(zf);
  free(zpath);
  char *nth = compute_nth_filename("x", ".y", 12);
  test_eq_str(nth, "x.y.12");
  free(nth);
  size_t n, calls = 0;
  char **paths = list_dir(dir, only_txt, &calls, &n, NULL);
  test_eq_size(n, 2);
  test_eq_size(calls, 4);   /* a, x.txt, y.dat, z.txt */
  base = path_basename(paths[0]);
  test_eq_str(base, "x.txt");
  free(base);
  for (size_t i = 0; i < n; i++)
    free(paths[i]);
  free(paths);
  paths = list_dir(dir, NULL, NULL, &n, NULL);
  test_eq_size(n, 4);   /* the hidden file is not listed */
  for (size_t i = 0; i < n; i++) {
    if (!is_directory(paths[i]))
      unlink(paths[i]);
    free(paths[i]);
  }
  free(paths);
  char *hidden = strcat_alloc(dir, "/.hidden.txt");
  unlink(hidden);
  free(hidden);
  rmdir(sub);
  char *a = strcat_alloc(dir, "/a");
  rmdir(a);
  free(a);
  rmdir(dir);
  free(sub);
  free(dir);
}

static void argparse_test(void) {
  long value;
  test_assert(argparse_long("42", &value));
  test_eq_long(value, 42L);
  test_assert(argparse_long("-3", &value));
  test_eq_long(value, -3L);
  test_assert(!argparse_long("3abc", &value));
  test_assert(!argparse_long("", &value));
  test_assert(!argparse_long("99999999999999999999", &value));
  /* the value must fit in the destination */
  int i;
  long big = 1L << 40, minus1 = -1;
  test_assert(ARGPARSE_STORE(&i, 42L));
  test_eq_int(i, 42);
  test_assert(!ARGPARSE_STORE(&i, big));
  unsigned int u;
  test_assert(!ARGPARSE_STORE(&u, minus1));
  size_t z;
  test_assert(ARGPARSE_STORE(&z, big));
  test_assert(!ARGPARSE_STORE(&z, minus1));
  test_eq_str(BOOLSTR(false), "F");
  /* PARSE_LOG and PARSE_LOG_N are statements (e.g. in an if-else without
     braces); with no logger, the LOG_PARAM messages are not output */
  Logger *parse_logger = NULL;
  if (i == 42)
    PARSE_LOG("Value: %d\n", i);
  else
    PARSE_LOG_N("other\n");
}

static void bits_test(void) {
  unsigned long long x = 0;
  SET_BIT(x, 40);
  SET_BIT(x, 0);
  test_assert(x == ((1ULL << 40) | 1));
  test_assert(GET_BIT(x, 40) == 1);
  test_assert(GET_BIT(x, 39) == 0);
  UNSET_BIT(x, 40);
  test_assert(x == 1);
  bloom_filter_test();
}

/* elements of the sorting test: lines of SORT_LINE_LEN characters */
#define SORT_LINE_LEN 6
#define SORT_N_LINES 1000

static int cmp_lines(const void *a, const void *b) {
  return memcmp(*(char *const *)a, *(char *const *)b, SORT_LINE_LEN);
}

static void sort_lines(char **lines, size_t n) {
  qsort(lines, n, sizeof(char *), cmp_lines);
}

static void output_lines(FILE *file, char **lines, size_t n,
    Logger *logger) {
  (void)logger;
  for (size_t i = 0; i < n; i++)
    fprintf(file, "%.*s\n", SORT_LINE_LEN, lines[i]);
}

static bool get_line(char **line, size_t *size, FILE *file, Logger *logger) {
  return xgetline(line, size, file, logger);
}

static void write_line(FILE *file, char *line, Logger *logger) {
  (void)logger;
  fputs(line, file);
}

/* the merger compares pointers to the elements (char **) */
static int cmp_line_ptrs(const void *a, const void *b) {
  return strcmp(*(char *const *)a, *(char *const *)b);
}

static void file_sorter_test(const char *tmpdir) {
  char *path = strcat_alloc(tmpdir, "/cbaia-tests-sort");
  FILE *f = fopen(path, "w");
  test_assert(f != NULL);
  for (size_t i = 0; i < SORT_N_LINES; i++)
    fprintf(f, "%06zu\n", (i * 7919) % 10007);
  fclose(f);
  /* small buffer: several parts */
  size_t nfiles = sort_file_parts(path, NULL, path, ".sorted", 1024, 16,
      read_lines_in_buffer, sort_lines, output_lines, NULL);
  test_assert(nfiles > 1);
  merge_n_files(nfiles, path, ".sorted", false, get_line, write_line,
                cmp_line_ptrs, NULL);
  char *sorted = strcat_alloc(path, ".sorted");
  f = fopen(sorted, "r");
  test_assert(f != NULL);
  char *line = NULL, prev[SORT_LINE_LEN + 2] = "";
  size_t size = 0, n = 0;
  while (xgetline(&line, &size, f, NULL)) {
    test_assert(strcmp(prev, line) <= 0);
    strcpy(prev, line);
    n++;
  }
  test_eq_size(n, SORT_N_LINES);
  fclose(f);
  unlink(sorted);
  /* a single part, to the standard output (redirected to a file): the part
     is copied, and closed (the next free descriptor is the same) */
  nfiles = sort_file_parts(path, NULL, path, ".sorted", 1 << 20, 16,
      read_lines_in_buffer, sort_lines, output_lines, NULL);
  test_eq_size(nfiles, 1);
  fflush(stdout);
  int saved_stdout = dup(1);
  f = fopen(sorted, "w");
  test_assert(f != NULL);
  dup2(fileno(f), 1);
  fclose(f);
  int fd_before = dup(0);
  close(fd_before);
  merge_n_files(nfiles, path, ".sorted", true, get_line, write_line,
                cmp_line_ptrs, NULL);
  int fd_after = dup(0);
  close(fd_after);
  fflush(stdout);
  dup2(saved_stdout, 1);
  close(saved_stdout);
  test_eq_int(fd_after, fd_before);
  f = fopen(sorted, "r");
  test_assert(f != NULL);
  n = 0;
  prev[0] = '\0';
  while (xgetline(&line, &size, f, NULL)) {
    test_assert(strcmp(prev, line) <= 0);
    strcpy(prev, line);
    n++;
  }
  test_eq_size(n, SORT_N_LINES);
  fclose(f);
  unlink(sorted);
  /* one line per part: more parts than are merged at once (512) */
  nfiles = sort_file_parts(path, NULL, path, ".sorted", SORT_LINE_LEN + 3,
      16, read_lines_in_buffer, sort_lines, output_lines, NULL);
  test_eq_size(nfiles, SORT_N_LINES);
  merge_n_files(nfiles, path, ".sorted", false, get_line, write_line,
                cmp_line_ptrs, NULL);
  f = fopen(sorted, "r");
  test_assert(f != NULL);
  n = 0;
  prev[0] = '\0';
  while (xgetline(&line, &size, f, NULL)) {
    test_assert(strcmp(prev, line) <= 0);
    strcpy(prev, line);
    n++;
  }
  test_eq_size(n, SORT_N_LINES);
  fclose(f);
  /* the temporary files were deleted */
  char *part = compute_nth_filename(path, ".sorted", 0);
  test_assert(!file_exists(part, NULL));
  free(part);
  part = compute_nth_filename(path, ".sorted", SORT_N_LINES);
  test_assert(!file_exists(part, NULL));
  free(part);
  free(line);
  unlink(sorted);
  unlink(path);
  free(sorted);
  free(path);
}

/* a test group: run or, for the tests which need the directory for
   temporary files, run_tmpdir */
typedef struct {
  char option;
  const char *description;
  void (*run)(void);
  void (*run_tmpdir)(const char *tmpdir);
  bool selected;
} test_group_t;

static test_group_t test_groups[] = {
  {'s', "string utilities and string buffer", strutils_test, NULL, false},
  {'N', "numeric utilities", numutils_test, NULL, false},
  {'k', "file and directory utilities", NULL, file_utils_test, false},
  {'a', "argument parsing", argparse_test, NULL, false},
  {'o', "progress display", progress_test, NULL, false},
  {'z', "gzip compressor", NULL, gzip_compressor_test, false},
  {'c', "BGZF compressor", NULL, bgzf_compressor_test, false},
  {'b', "bit operations and bloom filter", bits_test, NULL, false},
  {'e', "set", set_test, NULL, false},
  {'p', "heap", heap_test, NULL, false},
  {'y', "HyperLogLog", hll_test, NULL, false},
  {'f', "file sorting and merging", NULL, file_sorter_test, false},
  {'u', "sequence utilities", sequtil_test, NULL, false},
  {'E', "edit distance of short sequences", edit_distance_test, NULL, false},
  {'r', "FASTQ reader", NULL, fastq_reader_test, false},
  {'x', "FASTQ writer", fastq_writer_test, NULL, false},
  {'m', "BAM reader", NULL, bam_reader_test, false},
  {'n', "uBAM writer", NULL, ubam_writer_test, false},
  {'S', "SAM tags", sam_tags_test, NULL, false},
  {'A', "SAM header", NULL, sam_header_test, false},
  {'j', "MinKNOW read metadata", minknow_metadata_test, NULL, false},
  {'d', "reads reader (FASTQ or BAM)", NULL, reads_reader_test, false},
  {'t', "output files", NULL, outfiles_test, false},
  {'w', "barcode whitelist", NULL, whitelist_test, false},
  {'l', "cell caller", NULL, cell_caller_test, false},
  {'i', "library complexity", complexity_test, NULL, false},
  {'5', "MD5 digest", NULL, md5_test, false},
  {'W', "HTML writer", html_writer_test, NULL, false},
  {'M', "HTML templates", html_template_test, NULL, false},
  {'H', "HTML reports", html_page_test, NULL, false},
};
#define N_TEST_GROUPS (sizeof(test_groups) / sizeof(test_groups[0]))

static void print_purpose(FILE *stream, const char *progname) {
  fprintf(stream, "Unit tests of libcbaia\n"
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
        " -T, --tmpdir <dir> : directory for temporary files\n"
        "                      (default: " DEFAULT_TMPDIR ")\n"
        "\n"
        COMMON_OPTS_HELP, stream);
}

void parseargs(int argc, char *argv[],
    Logger **logger)
{
  PARSE_START;
  int c = 0;
  bool any = false;
  /* the other options, then one letter for each test group */
  char optstring[sizeof(COMMON_OPTS_STR "T:") + N_TEST_GROUPS] =
    COMMON_OPTS_STR "T:";
  size_t len = strlen(optstring);
  for (size_t i = 0; i < N_TEST_GROUPS; i++) {
    /* the test options must not be used by other options (checked also
       without assertions, as the tests are compiled with -DNDEBUG) */
    if (strchr(optstring, test_groups[i].option) != NULL)
      fatal_printf("Test option -%c used twice\n", test_groups[i].option);
    if (len + 1 >= sizeof(optstring))
      fatal_printf("Too many test options\n");
    optstring[len++] = test_groups[i].option;
    optstring[len] = '\0';
  }
  static struct option longopts[] = {
    {"tmpdir", required_argument, NULL, 'T'},
    COMMON_LONGOPTS
    LONGOPTS_END
  };
  opterr = 0;
  while ((c = getopt_long(argc, argv, optstring, longopts, NULL)) != -1) {
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
        opt_tmpdir = optarg;
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
  for (size_t i = 0; i < N_TEST_GROUPS; i++)
    if (test_groups[i].selected) {
      lprintf(logger, LOG_PHASE, "Testing %s...\n",
              test_groups[i].description);
      if (test_groups[i].run_tmpdir != NULL)
        test_groups[i].run_tmpdir(opt_tmpdir);
      else
        test_groups[i].run();
    }
  lprintf(logger, LOG_SUCCESS, "All tests passed\n");
  delete_logger(logger);
  return EXIT_SUCCESS;
}
