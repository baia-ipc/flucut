/* (c) 2026 Giorgio Gonnella */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "alloc.h"
#include "file.h"
#include "group-state.h"
#include "report.h"
#include "strutils.h"

/* the path of a state file of a group (to be freed by the caller) */
static char* state_path(const char *outdir, const char *group,
    const char *sfx) {
  char *dir = path_join(outdir, GROUP_STATE_DIR);
  char *name = strcat_alloc(group, sfx);
  char *path = path_join(dir, name);
  free(name);
  free(dir);
  return path;
}

void group_state_save(const char *outdir, const char *group,
    const counter_t *counts, char *const *files, size_t n_files,
    Logger *logger) {
  char *dir = path_join(outdir, GROUP_STATE_DIR);
  xmakedirs(dir, logger);
  free(dir);
  char *path = state_path(outdir, group, GROUP_COUNTS_SFX);
  FILE *f = xfopen(path, REPORT_TMP_SFX, "w", "state file", logger);
  counter_save(counts, f);
  xfclose(&f, logger);
  xfrename(path, REPORT_TMP_SFX, path, NULL, logger);
  free(path);
  path = state_path(outdir, group, GROUP_FILES_SFX);
  f = xfopen(path, REPORT_TMP_SFX, "w", "state file", logger);
  for (size_t i = 0; i < n_files; i++)
    fprintf(f, "%s\n", files[i]);
  xfclose(&f, logger);
  xfrename(path, REPORT_TMP_SFX, path, NULL, logger);
  free(path);
}

char** group_state_load_files(const char *outdir, const char *group,
    size_t *n_files) {
  *n_files = 0;
  char *path = state_path(outdir, group, GROUP_FILES_SFX);
  FILE *f = fopen(path, "r");
  free(path);
  if (f == NULL)
    return NULL;
  char **files = NULL, *line = NULL;
  size_t cap = 0, line_cap = 0;
  ssize_t len;
  while ((len = getline(&line, &line_cap, f)) > 0) {
    if (line[len - 1] == '\n')
      line[--len] = '\0';
    if (len == 0)
      continue;
    if (*n_files == cap) {
      cap = cap == 0 ? 16 : cap * 2;
      XREALLOC(files, sizeof(char *) * cap);
    }
    files[(*n_files)++] = xstrdup(line);
  }
  free(line);
  fclose(f);
  return files;
}

void group_state_free_files(char **files, size_t n_files) {
  for (size_t i = 0; i < n_files; i++)
    free(files[i]);
  free(files);
}

bool group_state_load_counts(const char *outdir, const char *group,
    counter_t *counts) {
  char *path = state_path(outdir, group, GROUP_COUNTS_SFX);
  FILE *f = fopen(path, "r");
  free(path);
  if (f == NULL)
    return false;
  bool ok = counter_load(counts, f);
  fclose(f);
  return ok;
}

void group_state_remove(const char *outdir, const char *group,
    Logger *logger) {
  const char *sfx[] = {GROUP_FILES_SFX, GROUP_COUNTS_SFX};
  for (size_t i = 0; i < 2; i++) {
    char *path = state_path(outdir, group, sfx[i]);
    LCHECK(unlink(path) == 0 || errno == ENOENT, logger,
           "Deleting file '%s' failed, %s\n", path, strerror(errno));
    free(path);
  }
}
