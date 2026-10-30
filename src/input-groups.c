/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "alloc.h"
#include "file.h"
#include "input-groups.h"
#include "reads-reader.h"
#include "strutils.h"
#include "test.h"

/* list_dir() filters: the barcode subdirectories, the files of reads */
static bool is_group_dir(const char *name, const char *path, void *data) {
  (void)data;
  return strncmp(name, INPUT_GROUPS_DIR_PFX,
                 strlen(INPUT_GROUPS_DIR_PFX)) == 0 && is_directory(path);
}

static bool is_reads_file(const char *name, const char *path, void *data) {
  (void)data;
  return reads_reader_is_reads_filename(name) && !is_directory(path);
}

static input_group_t *add_group(input_group_t *groups, size_t *n,
    char *name, char **files, size_t n_files) {
  XREALLOC(groups, sizeof(input_group_t) * (*n + 1));
  groups[*n].name = name;
  groups[*n].files = files;
  groups[*n].n_files = n_files;
  (*n)++;
  return groups;
}

/* the group of a file given directly */
static input_group_t *add_file(input_group_t *groups, size_t *n,
    const char *input, Logger *logger) {
  LCHECK(reads_reader_is_reads_filename(input), logger, "%s: not a file "
         "of reads (.fastq, .fq, optionally .gz, or .bam)\n", input);
  char *name = path_basename(input);
  name[strlen(name) - reads_reader_filename_suffix_len(name)] = '\0';
  char **files;
  XMALLOC(files, sizeof(char *));
  files[0] = xstrdup(input);
  return add_group(groups, n, name, files, 1);
}

input_group_t *input_groups_new(char *const *inputs, size_t n_inputs,
    size_t *n_groups, bool allow_empty, bool warn_ignored, Logger *logger) {
  assert(inputs != NULL || n_inputs == 0);
  assert(n_groups != NULL);
  input_group_t *groups = NULL;
  *n_groups = 0;
  for (size_t i = 0; i < n_inputs; i++) {
    const char *input = inputs[i];
    if (!is_directory(input)) {
      groups = add_file(groups, n_groups, input, logger);
      continue;
    }
    size_t n_subdirs, n_files;
    char **subdirs = list_dir(input, is_group_dir, NULL, &n_subdirs, logger);
    char **files = list_dir(input, is_reads_file, NULL, &n_files, logger);
    if (n_subdirs == 0)
      groups = add_group(groups, n_groups, path_basename(input), files,
                         n_files);
    else {
      /* the groups are the barcode directories: the files beside them
         are not read */
      if (n_files > 0 && warn_ignored)
        lprintf(logger, LOG_WARN, "The input directory %s contains "
                "barcode directories and %zu FASTQ or BAM file%s outside "
                "them (e.g. %s): only the barcode directories are "
                "processed\n", input, n_files, n_files > 1 ? "s" : "",
                files[0]);
      for (size_t j = 0; j < n_files; j++)
        free(files[j]);
      free(files);
    }
    for (size_t j = 0; j < n_subdirs; j++) {
      char **sub_files = list_dir(subdirs[j], is_reads_file, NULL,
                                  &n_files, logger);
      groups = add_group(groups, n_groups, path_basename(subdirs[j]),
                         sub_files, n_files);
      free(subdirs[j]);
    }
    free(subdirs);
  }
  for (size_t i = 0; i < *n_groups; i++) {
    if (groups[i].n_files == 0 && !allow_empty)
      lprintf(logger, LOG_WARN, "No FASTQ or BAM files in the input group "
              "%s\n", groups[i].name);
    for (size_t j = 0; j < i; j++)
      LCHECK(strcmp(groups[i].name, groups[j].name) != 0, logger,
             "The input group name %s is not unique (both groups would be "
             "written to %s" INPUT_GROUPS_OUTPUT_SFX ")\n", groups[i].name,
             groups[i].name);
  }
  return groups;
}

void input_groups_delete(input_group_t *groups, size_t n_groups) {
  for (size_t i = 0; i < n_groups; i++) {
    for (size_t j = 0; j < groups[i].n_files; j++)
      free(groups[i].files[j]);
    free(groups[i].files);
    free(groups[i].name);
  }
  free(groups);
}

char *input_group_output_path(const char *outdir, const char *group) {
  assert(outdir != NULL);
  assert(group != NULL);
  char *filename = strcat_alloc(group, INPUT_GROUPS_OUTPUT_SFX);
  char *path = path_join(outdir, filename);
  free(filename);
  return path;
}

/* --- unit tests --- */

/* the paths created by the test, in the order of their creation */
#define TEST_MAX_PATHS 16

typedef struct {
  char *paths[TEST_MAX_PATHS];
  size_t n;
} test_tree_t;

/* create <tmpdir>/<name>: a directory if <name> ends with a slash,
   otherwise an empty file */
static void test_create(test_tree_t *t, const char *tmpdir,
    const char *name) {
  assert(t->n < TEST_MAX_PATHS);
  char *path = path_join(tmpdir, name);
  size_t len = strlen(path);
  if (path[len - 1] == '/') {
    path[len - 1] = '\0';
    test_assert(mkdir(path, 0700) == 0);
  } else {
    FILE *f = fopen(path, "w");
    test_assert(f != NULL);
    fclose(f);
  }
  t->paths[t->n++] = path;
}

static void test_group(const input_group_t *g, const char *name,
    size_t n_files, const char *tmpdir, const char *first) {
  test_eq_str(g->name, name);
  test_eq_size(g->n_files, n_files);
  if (first != NULL) {
    char *path = path_join(tmpdir, first);
    test_eq_str(g->files[0], path);
    free(path);
  }
}

void input_groups_test(const char *tmpdir) {
  test_tree_t t = {.n = 0};
  const char *tree[] = {"ig-s1.fastq.gz", "ig-run/", "ig-run/fastq_pass/",
      "ig-run/fastq_pass/barcode01/", "ig-run/fastq_pass/barcode01/a.fq",
      "ig-run/fastq_pass/barcode02/", "ig-run/fastq_pass/barcode02/c.bam",
      "ig-run/fastq_pass/barcode02/b.fastq.gz",
      "ig-run/fastq_pass/barcode02/notes.txt",
      "ig-run/fastq_pass/unclassified/",
      "ig-run/fastq_pass/unclassified/x.fastq",
      "ig-run/fastq_pass/beside.fastq", "ig-run/plain/",
      "ig-run/plain/d.fastq", "ig-run/empty/", NULL};
  for (const char **n = tree; *n != NULL; n++)
    test_create(&t, tmpdir, *n);
  char *s1 = path_join(tmpdir, "ig-s1.fastq.gz"),
       *pass = path_join(tmpdir, "ig-run/fastq_pass/"),
       *plain = path_join(tmpdir, "ig-run/plain");
  char *inputs[] = {s1, pass, plain};
  size_t n;
  input_group_t *g = input_groups_new(inputs, 3, &n, false, false, NULL);
  /* the file given directly, the barcode directories (the files beside
     them and the other subdirectories are not read), the directory */
  test_eq_size(n, 4);
  test_group(g, "ig-s1", 1, tmpdir, "ig-s1.fastq.gz");
  test_group(g + 1, "barcode01", 1, tmpdir,
             "ig-run/fastq_pass/barcode01/a.fq");
  test_group(g + 2, "barcode02", 2, tmpdir,
             "ig-run/fastq_pass/barcode02/b.fastq.gz");
  test_group(g + 3, "plain", 1, tmpdir, "ig-run/plain/d.fastq");
  input_groups_delete(g, n);
  /* a directory without files of reads: an empty group */
  char *empty = path_join(tmpdir, "ig-run/empty/");
  g = input_groups_new(&empty, 1, &n, true, true, NULL);
  test_eq_size(n, 1);
  test_group(g, "empty", 0, tmpdir, NULL);
  input_groups_delete(g, n);
  char *out = input_group_output_path("out/", "barcode01");
  test_eq_str(out, "out/barcode01.fastq.gz");
  free(out);
  free(empty);
  free(s1);
  free(pass);
  free(plain);
  while (t.n > 0) {
    t.n--;
    test_assert(remove(t.paths[t.n]) == 0);
    free(t.paths[t.n]);
  }
}
