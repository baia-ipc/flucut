/* (c) 2026 Giorgio Gonnella */

/*
   Writing the JSON summary: nested objects and arrays (blocks) with one
   member on each line, indented by the nesting level; the values of the
   members are written on their line by the caller (lists on one line:
   see json_sep()). The writer keeps, for each open block, whether it has
   a member yet, for the commas between the members.
*/

#ifndef JSON_WRITER_H
#define JSON_WRITER_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

#define JSON_WRITER_MAX_NESTING 8

typedef struct {
  FILE *f;
  int level;  /* number of open blocks */
  bool empty[JSON_WRITER_MAX_NESTING];
} json_writer_t;

/* Open a block: <c> is { or [ */
void json_open_block(json_writer_t *w, char c);

/* Start a member of the open block: on a new line (after a comma, unless
   the first); with "<k>": (k != NULL: object) */
void json_member(json_writer_t *w, const char *k);

/* Close the open block: <c> is } or ]; on a new line, or, if empty and
   <compact>, as {} or [] */
void json_close_block(json_writer_t *w, char c, bool compact);

/* The separator of the elements of a list on one line (none before the
   first: <first>, then set to false) */
void json_sep(FILE *f, bool *first);

/* A member with a count; with a share (a fraction of 1) */
void json_count_member(json_writer_t *w, const char *k, size_t n);
void json_share_member(json_writer_t *w, const char *k, double share);

#endif
