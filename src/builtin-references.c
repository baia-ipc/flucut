/* (c) 2026 Giorgio Gonnella */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "file.h"
#include "builtin-references.h"

/*
   The accessions of the built-in set: complete segments of influenza A
   from NCBI, the HA of every subtype (H1 to H16), the NA of every
   subtype (N1 to N9) and the other six segments of isolates of several
   subtypes. Each one begins and ends with the conserved termini of the
   segments (the targets of the Uni12 and Uni13 primers), so that the
   ends of the reads at a segment terminus are recognised as such, and
   each carries its subtype in its description. The script
   find_builtin_references.py in scripts/ searches NCBI for records
   which fulfil this and writes the list below. The subtype of a
   sample is the subtype of the HA and of the NA its reads are closest
   to, so every subtype needs its own reference; the other segments
   differ little between subtypes, and several isolates cover their
   diversity. Updating this list requires a
   new BUILTIN_REFERENCES_VERSION (the cache file is named after it, so
   an old cache is not reused).
*/
#define BUILTIN_REFERENCES_VERSION "v2"
static const char *const BUILTIN_ACCESSIONS =
  /* the HA of each subtype, H1 to H16 */
  "QB342831.1,QB057580.1,QB367364.1,QB139075.1,"
  "PV571614.1,QB057598.1,OL782158.1,PZ231346.1,"
  "QB023536.1,QB024178.1,QB023748.1,QB024967.1,"
  "PV450481.1,OP144199.1,QB024617.1,QB024068.1,"
  /* the NA of each subtype, N1 to N9 */
  "QB024910.1,QB367366.1,QB023727.1,PZ229867.1,"
  "QB024743.1,QB024853.1,QB024172.1,QB025125.1,"
  "QB057582.1,"
  /* the other segments (PB2, PB1, PA, NP, MP, NS) of isolates of
     different subtypes, years and hosts */
  /* H1N1 */
  "QB367363.1,QB367331.1,QB342828.1,QB367360.1,"
  "QB367357.1,QB367329.1,"
  /* H3N2 */
  "QB367355.1,QB367354.1,PZ944660.1,QB367352.1,"
  "QB367365.1,QB367367.1,"
  /* H5N1 */
  "PQ064247.4,PZ598243.1,PZ598252.1,PZ599870.1,"
  "QB332476.1,QB024702.1,"
  /* H5N8 */
  "PV410901.1,EU743133.2,PV410903.1,EU743131.2,"
  "PZ598696.1,PZ599841.1,"
  /* H7N9 */
  "PQ756676.1,PQ756677.1,PQ756673.1,PQ756675.1,"
  "PQ756678.1,PQ756679.1,"
  /* H9N2 */
  "QB023906.1,PZ768253.1,ON328346.1,PZ768255.1,"
  "QB214291.1,QB215690.1";

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
