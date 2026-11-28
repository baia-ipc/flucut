/* (c) 2026 Giorgio Gonnella */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "test.h"
#include "user-config.h"

#define CONFIG_LINE_SIZE 4096

/* the path of the configuration file (a static buffer, valid until the
   next call); NULL if neither XDG_CONFIG_HOME nor HOME is set */
static const char* config_path(const char *program) {
  static char path[4096];
  const char *xdg = getenv("XDG_CONFIG_HOME");
  if (xdg != NULL && xdg[0] != '\0') {
    snprintf(path, sizeof(path), "%s/%s/config", xdg, program);
    return path;
  }
  const char *home = getenv("HOME");
  if (home != NULL && home[0] != '\0') {
    snprintf(path, sizeof(path), "%s/.config/%s/config", home, program);
    return path;
  }
  return NULL;
}

/* the text without the spaces and tabs at its ends (in place) */
static char* trim(char *s) {
  while (*s == ' ' || *s == '\t')
    s++;
  char *end = s + strlen(s);
  while (end > s && (end[-1] == ' ' || end[-1] == '\t' || end[-1] == '\n' ||
                     end[-1] == '\r'))
    end--;
  *end = '\0';
  return s;
}

/* the value of <key> in the file <path>, written to <value> (size
   <size>); false if the file cannot be read or the key is not there */
static bool value_in_file(const char *path, const char *key, char *value,
    size_t size) {
  FILE *f = fopen(path, "r");
  if (f == NULL)
    return false;
  char line[CONFIG_LINE_SIZE];
  bool found = false;
  while (!found && fgets(line, sizeof(line), f) != NULL) {
    char *hash = strchr(line, '#');
    if (hash != NULL)
      *hash = '\0';
    char *equals = strchr(line, '=');
    if (equals == NULL)
      continue;
    *equals = '\0';
    if (strcmp(trim(line), key) != 0)
      continue;
    char *text = trim(equals + 1);
    if (text[0] != '\0') {
      snprintf(value, size, "%s", text);
      found = true;
    }
  }
  fclose(f);
  return found;
}

const char* user_config_value(const char *program, const char *key,
    Logger *logger) {
  static char value[CONFIG_LINE_SIZE];
  const char *path = config_path(program);
  if (path == NULL)
    return NULL;
  if (!value_in_file(path, key, value, sizeof(value)))
    return NULL;
  lprintf(logger, LOG_DEBUG, "%s of the configuration file %s: %s\n", key,
          path, value);
  return value;
}

void user_config_test(void) {
  const char *name = "fc-test-config";
  FILE *f = fopen(name, "w");
  test_assert(f != NULL);
  fputs("# a comment\n"
        "  references = /tmp/refs.fasta   # where they are\n"
        "irma=/opt/IRMA\n"
        "empty =\n"
        "no equals sign\n", f);
  fclose(f);
  char value[CONFIG_LINE_SIZE];
  test_assert(value_in_file(name, "references", value, sizeof(value)));
  test_eq_str(value, "/tmp/refs.fasta");
  test_assert(value_in_file(name, "irma", value, sizeof(value)));
  test_eq_str(value, "/opt/IRMA");
  test_assert(!value_in_file(name, "empty", value, sizeof(value)));
  test_assert(!value_in_file(name, "other", value, sizeof(value)));
  test_assert(!value_in_file("fc-test-no-such-file", "references", value,
                           sizeof(value)));
  remove(name);
}
