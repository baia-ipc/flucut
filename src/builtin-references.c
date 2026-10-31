/* (c) 2026 Giorgio Gonnella */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "file.h"
#include "builtin-references.h"

/*
   The accessions of the built-in set: one complete segment (with the
   UTRs) per segment (PB2, PB1, PA, HA, NP, NA, MP, NS) of each subtype,
   from NCBI; each is verified to carry its subtype in its description.
   Updating this list requires a new BUILTIN_REFERENCES_VERSION (the cache
   file is named after it, so an old cache is not reused).
*/
#define BUILTIN_REFERENCES_VERSION "v1"
static const char *const BUILTIN_ACCESSIONS =
  /* H1N1 (pdm09) */
  "QB169385.1,QB169376.1,QB169383.1,QB169386.1,"
  "QB169381.1,QB169363.1,QB230774.1,QB169382.1,"
  /* H3N2 */
  "QB100355.1,QB100356.1,QB100316.1,QB100352.1,"
  "QB100304.1,QB100340.1,QB100350.1,QB211866.1,"
  /* H5N1 */
  "QB222556.1,QB222557.1,QB222558.1,QB222559.1,"
  "QB222560.1,QB222561.1,QB222562.1,QB222563.1,"
  /* H7N9 */
  "PQ756676.1,PQ756677.1,PQ756673.1,PQ756680.1,"
  "PQ756675.1,OQ820951.1,PQ756678.1,PQ756679.1,"
  /* H9N2 */
  "QB023906.1,QB023903.1,QB023392.1,QB001843.1,"
  "QB023389.1,QB024010.1,QB214291.1,QB215583.1";

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

/* download the accessions to <path> (via curl, which honors HTTPS_PROXY);
   false on failure */
static bool download(const char *path, Logger *logger) {
  size_t n = strlen(path) + 16;
  char *tmp_path = malloc(n);
  if (tmp_path == NULL)
    return false;
  snprintf(tmp_path, n, "%s.part", path);
  char *cmd;
  size_t cmd_size = strlen(BUILTIN_ACCESSIONS) + strlen(tmp_path) + 512;
  cmd = malloc(cmd_size);
  if (cmd == NULL) {
    free(tmp_path);
    return false;
  }
  snprintf(cmd, cmd_size,
      "curl -sSf --retry 3 --retry-delay 2 -o '%s' "
      "'%s?db=nuccore&id=%s&rettype=fasta&retmode=text'",
      tmp_path, NCBI_EFETCH, BUILTIN_ACCESSIONS);
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

const char* builtin_references_path(Logger *logger) {
  const char *override = getenv("FLUCUT_BUILTIN_REFERENCES");
  if (override != NULL && override[0] != '\0')
    return override;
  static char path[4096];
  const char *dir = cache_dir();
  if (dir == NULL) {
    lprintf(logger, LOG_WARN, "no cache directory (set XDG_CACHE_HOME or "
            "HOME) for the built-in references; running without references "
            "(use -r <file>, or -r off)\n");
    return NULL;
  }
  snprintf(path, sizeof(path), "%s/flu-a-references-%s.fasta", dir,
           BUILTIN_REFERENCES_VERSION);
  if (file_ok(path))
    return path;
  xmakedirs(dir, logger);
  if (download(path, logger))
    return path;
  lprintf(logger, LOG_WARN, "could not download the built-in references "
          "(no network, or curl missing); running without references "
          "(use -r <file>, or -r off to silence this)\n");
  return NULL;
}
