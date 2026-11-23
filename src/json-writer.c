/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include "json-writer.h"
#include "strutils.h"

static void indent(FILE *f, int level) {
  for (int i = 0; i < level; i++)
    fputs("  ", f);
}

/* open a block: <c> is { or [ */
void json_open_block(json_writer_t *w, char c) {
  assert(w->level < JSON_WRITER_MAX_NESTING);
  fputc(c, w->f);
  w->empty[w->level++] = true;
}

/* start a member of the open block: on a new line (after a comma, unless
   the first); with "<k>": (k != NULL: object) */
void json_member(json_writer_t *w, const char *k) {
  assert(w->level > 0);
  fputs(w->empty[w->level - 1] ? "\n" : ",\n", w->f);
  w->empty[w->level - 1] = false;
  indent(w->f, w->level);
  if (k != NULL) {
    fputs_json_string(k, w->f);
    fputs(": ", w->f);
  }
}

/* close the open block: <c> is } or ]; on a new line, or, if empty and
   <compact>, as {} or [] */
void json_close_block(json_writer_t *w, char c, bool compact) {
  assert(w->level > 0);
  w->level--;
  if (!w->empty[w->level] || !compact) {
    fputc('\n', w->f);
    indent(w->f, w->level);
  }
  fputc(c, w->f);
}

/* the separator of the elements of a list on one line (none before the
   first: <first>, then set to false) */
void json_sep(FILE *f, bool *first) {
  if (!*first)
    fputs(", ", f);
  *first = false;
}

/* a member with a count */
void json_count_member(json_writer_t *w, const char *k, size_t n) {
  json_member(w, k);
  fprintf(w->f, "%zu", n);
}

/* a share (a fraction of 1) */
void json_share_member(json_writer_t *w, const char *k, double share) {
  json_member(w, k);
  fprintf(w->f, "%.4f", share);
}
