/* (c) 2017 Giorgio Gonnella */

#ifndef STRUTILS_H
#define STRUTILS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

/*
  Concatenate the two 0-term strings s1 and s2 in
  a newly allocated string.

  The user must remember to free the buffer containing the returned value.
*/
char *strcat_alloc(const char *s1, const char *s2);

/* Copy of the string s (allocated; the program is exited if the
   allocation fails, unlike strdup()) */
char *xstrdup(const char *s);

/* Whether the string s ends with sfx (and is longer than it) */
bool str_has_suffix(const char *s, const char *sfx);

/* Write s to f as a JSON string (with the quotes; the characters which
   must be escaped are escaped) */
void fputs_json_string(const char *s, FILE *f);

/* Write s to f, escaping the HTML special characters (& < > ") */
void fputs_html_escaped(const char *s, FILE *f);

/* Format a duration in seconds as [D d ]HH:MM:SS (rounded to the
   second; negative durations as 0) into buf, of the given size */
void format_duration(double seconds, char *buf, size_t size);

/* Parse a duration: a non-negative number with an optional unit s, m or h
   (e.g. "30s", "1.5h"; without unit: in units of default_unit seconds,
   e.g. 60 for minutes) into seconds (rounded); false if invalid */
bool parse_duration(const char *s, double default_unit, long *seconds);

/* Whether c is whitespace, as for the split() of Python (ASCII: space,
   \t, \n, \v, \f, \r and the separators \x1c-\x1f) */
static inline bool is_space(char c) {
  return c == ' ' || (c >= '\t' && c <= '\r') ||
         (c >= '\x1c' && c <= '\x1f');
}

#endif
