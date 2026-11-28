/* (c) 2026 Giorgio Gonnella */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "file.h"
#include "builtin-references.h"
#include "user-config.h"

/*
   The accessions of the built-in sets and the version of each, written
   by the build from metadata/references.yaml
   (scripts/references-table.awk). A reference is added or changed
   there; the version of a set must change with its list, as the cache
   file is named after it and an old cache would otherwise be reused.
*/
#include "references-table.h"

#define NCBI_EFETCH \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

/* the directory for the cache: $XDG_CACHE_HOME/flucut, else ~/.cache/flucut
   (a static buffer, valid until the next call); NULL if neither is set */
static const char* cache_dir(void) {
  static char dir[4096];
  const char *xdg = getenv("XDG_CACHE_HOME");
  if (xdg != NULL && xdg[0] != '\0') {
    snprintf(dir, sizeof(dir), "%s/flucut", xdg);
    return dir;
  }
  const char *home = getenv("HOME");
  if (home != NULL && home[0] != '\0') {
    snprintf(dir, sizeof(dir), "%s/.cache/flucut", home);
    return dir;
  }
  return NULL;
}

/* whether <path> is an existing non-empty file */
static bool file_ok(const char *path) {
  struct stat st;
  return stat(path, &st) == 0 && S_ISREG(st.st_mode) && st.st_size > 0;
}

/* download the accessions <accessions> to <path> (via curl, which
   honors HTTPS_PROXY); false on failure */
static bool download(const char *path, const char *accessions,
    Logger *logger) {
  size_t n = strlen(path) + 16;
  char *tmp_path = malloc(n);
  if (tmp_path == NULL)
    return false;
  snprintf(tmp_path, n, "%s.part", path);
  char *cmd;
  size_t cmd_size = strlen(accessions) + strlen(tmp_path) + 512;
  cmd = malloc(cmd_size);
  if (cmd == NULL) {
    free(tmp_path);
    return false;
  }
  snprintf(cmd, cmd_size,
      "curl -sSf --retry 3 --retry-delay 2 -o '%s' "
      "'%s?db=nuccore&id=%s&rettype=fasta&retmode=text'",
      tmp_path, NCBI_EFETCH, accessions);
  lprintf(logger, LOG_PHASE,
          "Downloading the built-in reference set (first use)...\n");
  int status = system(cmd);
  free(cmd);
  bool ok = status == 0 && file_ok(tmp_path);
  if (ok && rename(tmp_path, path) != 0)
    ok = false;
  if (!ok)
    remove(tmp_path);
  free(tmp_path);
  return ok;
}

/*
   The accessions of the set of the types <set> (both types: those of
   influenza A and then those of influenza B, in one download), and the
   name of its cache file; the accessions are returned and <name> is
   filled. The caller frees the result.
*/
static char* set_accessions(virus_set_t set, char *name, size_t size) {
  bool a = virus_set_has(set, VIRUS_A), b = virus_set_has(set, VIRUS_B);
  size_t n = (a ? strlen(BUILTIN_ACCESSIONS_A) : 0) +
             (b ? strlen(BUILTIN_ACCESSIONS_B) : 0) + 2;
  char *accessions = malloc(n);
  if (accessions == NULL)
    return NULL;
  snprintf(accessions, n, "%s%s%s", a ? BUILTIN_ACCESSIONS_A : "",
           a && b ? "," : "", b ? BUILTIN_ACCESSIONS_B : "");
  if (a && b)
    snprintf(name, size, "flu-ab-references-%s-%s.fasta",
             BUILTIN_REFERENCES_A_VERSION, BUILTIN_REFERENCES_B_VERSION);
  else if (b)
    snprintf(name, size, "flu-b-references-%s.fasta",
             BUILTIN_REFERENCES_B_VERSION);
  else
    snprintf(name, size, "flu-a-references-%s.fasta",
             BUILTIN_REFERENCES_A_VERSION);
  return accessions;
}

const char* builtin_references_path(virus_set_t set, Logger *logger) {
  const char *override = getenv("FLUCUT_BUILTIN_REFERENCES");
  if (override != NULL && override[0] != '\0')
    return override;
  static char path[4096];
  /* a set of the configuration file, e.g. one built from the seeds of
     IRMA by scripts/irma_references.py (see user-config.h) */
  const char *configured = user_config_value("flucut", "references",
                                              logger);
  if (configured != NULL) {
    snprintf(path, sizeof(path), "%s", configured);
    if (file_ok(path)) {
      lprintf(logger, LOG_INFO, "references of the configuration file: "
              "%s\n", path);
      return path;
    }
    lprintf(logger, LOG_WARN, "the references of the configuration file "
            "(%s) are not there; using the built-in set\n", path);
  }
  const char *dir = cache_dir();
  if (dir == NULL) {
    lprintf(logger, LOG_WARN, "no cache directory (set XDG_CACHE_HOME or "
            "HOME) for the built-in references; running without references "
            "(use -r <file>, or -r off)\n");
    return NULL;
  }
  char name[64];
  char *accessions = set_accessions(set, name, sizeof(name));
  if (accessions == NULL)
    return NULL;
  snprintf(path, sizeof(path), "%s/%s", dir, name);
  if (file_ok(path)) {
    free(accessions);
    return path;
  }
  xmakedirs(dir, logger);
  bool ok = download(path, accessions, logger);
  free(accessions);
  if (ok)
    return path;
  lprintf(logger, LOG_WARN, "could not download the built-in references "
          "(no network, or curl missing); running without references "
          "(use -r <file>, or -r off to silence this)\n");
  return NULL;
}
