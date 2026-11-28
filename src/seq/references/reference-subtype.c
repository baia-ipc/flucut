/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "alloc.h"
#include "reference-subtype.h"
#include "strutils.h"
#include "virus.h"

/*
   The lineage of influenza B in a FASTA description: "Victoria" or
   "Yamagata" where the description says so ("Victoria lineage",
   "lineage Yamagata"), and "B" where it does not, so that the
   references whose lineage is not known are one group. Never NULL.

   The word alone is not enough: Victoria and Yamagata are also the
   places of the strain names B/Victoria/2/87 and B/Yamagata/16/88,
   which name the lineages after them, and a virus of either lineage
   isolated in Victoria or in Yamagata carries that place in its name.
   Most records do not say which lineage they are, so most references
   of influenza B are the group "B": the lineage of a sample is better
   called downstream, from its assembled HA.
*/
static char* parse_lineage(const char *desc) {
  static const char *const names[] = {"Victoria", "Yamagata"};
  static const char key[] = "lineage";
  const size_t key_len = sizeof(key) - 1;
  for (size_t i = 0; i < 2 && desc != NULL; i++) {
    size_t len = strlen(names[i]);
    for (const char *p = desc; (p = strcasestr(p, names[i])) != NULL;
         p += len) {
      /* "<name> lineage" or "lineage <name>" (any separators between) */
      const char *after = p + len;
      while (*after == ' ' || *after == '-' || *after == '_')
        after++;
      if (strncasecmp(after, key, key_len) == 0)
        return xstrdup(names[i]);
      if (p - desc >= (long)key_len + 1) {
        const char *before = p - 1;
        while (before > desc &&
               (*before == ' ' || *before == '-' || *before == '_'))
          before--;
        if (before - desc >= (long)key_len - 1 &&
            strncasecmp(before - key_len + 1, key, key_len) == 0)
          return xstrdup(names[i]);
      }
    }
  }
  return xstrdup("B");
}

char* reference_parse_subtype(const char *desc) {
  if (desc == NULL)
    return NULL;
  if (virus_description_is_b(desc))
    return parse_lineage(desc);
  for (const char *p = desc; (p = strchr(p, 'H')) != NULL; p++) {
    const char *q = p + 1;
    if (!isdigit((unsigned char)*q))
      continue;
    while (isdigit((unsigned char)*q))
      q++;
    if (*q != 'N')
      continue;
    const char *r = q + 1;
    if (!isdigit((unsigned char)*r))
      continue;
    while (isdigit((unsigned char)*r))
      r++;
    size_t n = (size_t)(r - p);
    char *s;
    XMALLOC(s, n + 1);
    for (size_t i = 0; i < n; i++)
      s[i] = (char)toupper((unsigned char)p[i]);
    s[n] = '\0';
    return s;
  }
  return NULL;
}

bool reference_description_subtype(const char *description, char *buf,
    size_t size) {
  assert(buf != NULL && size > 0);
  buf[0] = '\0';
  char *name = reference_parse_subtype(description);
  if (name == NULL)
    return false;
  snprintf(buf, size, "%s", name);
  free(name);
  return true;
}
