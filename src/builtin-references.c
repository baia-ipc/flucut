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
   to, so every subtype needs references of its own: several of them,
   as the lineages of one subtype can lie decades apart and one record
   does not stand for all of them (the votes of the references of a
   subtype are counted together, see auto-settings.h). The other
   segments differ little between subtypes, and several isolates cover
   their diversity. Updating this list requires a
   new BUILTIN_REFERENCES_VERSION (the cache file is named after it, so
   an old cache is not reused).
*/
#define BUILTIN_REFERENCES_VERSION "v3"
static const char *const BUILTIN_ACCESSIONS_A =
  /* the HA of each subtype, H1 to H16: four references of each,
     of different decades, hosts and places, as the lineages of a
     subtype can lie decades apart */
  /* H1 */ "QB342831.1,AF091316.1,KC894804.1,MF872838.1,"
  /* H2 */ "QB057580.1,KM213380.1,MG280380.1,KM885170.1,"
  /* H3 */ "QB367364.1,KP202378.1,MF872855.1,KP202377.1,"
  /* H4 */ "QB139075.1,GU066568.1,KX530513.2,GU066567.1,"
  /* H5 */ "PV571614.1,MN037415.1,MH276986.1,MW132374.1,"
  /* H6 */ "QB057598.1,MF575311.1,MW110064.1,MW110056.1,"
  /* H7 */ "OL782158.1,AY338460.1,MF978398.1,AY338459.1,"
  /* H8 */ "PZ231346.1,MG280398.1,MK995833.1,PZ229869.1,"
  /* H9 */ "QB023536.1,KP100808.1,MN037524.1,KC821026.2,"
  /* H10 */ "QB024178.1,GQ404728.2,KT887402.1,MG280186.1,"
  /* H11 */ "QB023748.1,JN798212.2,MT818622.1,KF542875.1,"
  /* H12 */ "QB024967.1,JX566034.1,MT818624.1,MT407041.1,"
  /* H13 */ "PV450481.1,KP033522.1,MF575207.1,KJ372720.1,"
  /* H14 */ "OP144199.1,KJ195668.1,MK928228.1,JN696316.2,"
  /* H15 */ "QB024617.1,MF147992.1,MF147328.1,MF146641.1,"
  /* H16 */ "QB024068.1,MK192316.1,MH158642.1,MZ015712.1,"
  /* the NA of each subtype, N1 to N9, likewise */
  /* N1 */ "QB024910.1,MF872852.1,QB342833.1,MW132361.1,"
  /* N2 */ "QB367366.1,KC894807.1,MN037526.1,KC345623.1,"
  /* N3 */ "QB023727.1,AY611526.1,KF885677.1,MW081048.1,"
  /* N4 */ "PZ229867.1,MG280500.1,MG280039.1,PP829796.1,"
  /* N5 */ "QB024743.1,M24740.1,MW132387.1,MT818788.1,"
  /* N6 */ "QB024853.1,MW186789.1,MF575219.1,QB023868.1,"
  /* N7 */ "QB024172.1,MW325955.1,MT406971.1,MK913511.1,"
  /* N8 */ "QB025125.1,KP202376.1,MW132362.1,KP202374.1,"
  /* N9 */ "QB057582.1,KF824506.1,MT773562.1,KF824505.1,"
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

/*
   The accessions of the built-in set of influenza B, found the same
   way (scripts/find_builtin_references.py --virus b). Influenza B has
   no subtypes: its diversity is two lineages, Victoria and Yamagata,
   which the records rarely name, so the set holds four references of
   each segment instead, of different places and decades, and the
   reference of each segment closest to the reads of a sample is chosen
   from them as for influenza A. Every record begins and ends with the
   conserved termini of influenza B, which are shorter than those of
   influenza A and continue into bases which differ between the
   segments (see primers.h): a reference without them would turn the
   ends of the reads at a segment terminus into internal priming.
   Updating this list requires a new BUILTIN_REFERENCES_B_VERSION.
*/
#define BUILTIN_REFERENCES_B_VERSION "v1"
static const char *const BUILTIN_ACCESSIONS_B =
  /* PB2 (segment 2 of influenza B) */
  "AF005737.1,AY504607.1,MN700020.1,PX491242.1,"
  /* PB1 (its segment 1) */
  "AY504608.1,MN700019.1,PZ093196.1,JX513083.1,"
  /* PA */
  "AF005738.1,KC866602.1,MN700021.1,PV226097.1,"
  /* HA */
  "EF626636.1,KX058884.1,PZ093207.1,PX494198.1,"
  /* NP */
  "MN700018.1,PX491246.1,EF626639.1,KC866605.1,"
  /* NA */
  "FJ766841.1,PZ093225.1,JQ396216.1,PX363335.1,"
  /* MP */
  "KC866607.1,MT040797.1,QB377036.1,KX351448.1,"
  /* NS */
  "KC866608.1,KM654633.1,PX363337.1,KC866606.1";

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
             BUILTIN_REFERENCES_VERSION, BUILTIN_REFERENCES_B_VERSION);
  else if (b)
    snprintf(name, size, "flu-b-references-%s.fasta",
             BUILTIN_REFERENCES_B_VERSION);
  else
    snprintf(name, size, "flu-a-references-%s.fasta",
             BUILTIN_REFERENCES_VERSION);
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
