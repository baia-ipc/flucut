/* (c) 2026 Giorgio Gonnella */

/*
   The pieces of markup shared by the parts of the HTML report: elements
   with a class, table cells, tiles, status badges, collapsed
   explanations, and the stacked bar, which shows how a total divides
   into a few parts (e.g. the reads by primers found). Their style is in
   html-report-style.html.
*/

#ifndef HTML_REPORT_MARKUP_H
#define HTML_REPORT_MARKUP_H

#include <stddef.h>
#include "html-writer.h"

/* Start an element of the class <cls> */
void markup_open(html_writer_t *w, const char *tag, const char *cls);

/* An element of the class <cls> with the text formatted as by printf */
void markup_textf(html_writer_t *w, const char *tag, const char *cls,
    const char *fmt, ...) __attribute__((format(printf, 4, 5)));

/* A small, muted paragraph, with the text formatted as by printf */
void markup_note(html_writer_t *w, const char *fmt, ...)
  __attribute__((format(printf, 2, 3)));

/* Start a table row (written on one line) */
void markup_open_row(html_writer_t *w);

/* A header row: the cells contain the given markup (constant strings), up
   to a NULL */
void markup_header_row(html_writer_t *w, ...);

/* a table cell with the text formatted as by printf */
#define markup_td(W, ...) html_tag_textf((W), "td", __VA_ARGS__)

/* A table cell without a value (a dash) */
void markup_td_none(html_writer_t *w);

/* A tile: a value (formatted as by printf) with its <label> and, under
   it, a smaller remark <sub> (NULL: none) */
void markup_tile(html_writer_t *w, const char *label, const char *sub,
    const char *fmt, ...) __attribute__((format(printf, 4, 5)));

/* the states shown by a badge */
typedef enum { MARKUP_GOOD, MARKUP_WARNING, MARKUP_CRITICAL,
               MARKUP_NEUTRAL } markup_state_t;

/* A badge: a state, as a symbol and a color, with a short text formatted
   as by printf ("": the symbol alone) */
void markup_badge(html_writer_t *w, markup_state_t state, const char *fmt,
    ...) __attribute__((format(printf, 3, 4)));

/* Start a part shown only when the reader opens it, under the line
   <summary>; end it with markup_details_end() */
void markup_details_start(html_writer_t *w, const char *summary);
void markup_details_end(html_writer_t *w);

/* a part of a stacked bar: <n> items, drawn in the color of the class
   <cls> */
typedef struct {
  const char *cls, *label;
  size_t n;
} markup_part_t;

/* A stacked bar titled <title>: the <n> <parts> as shares of their sum,
   which is given after the title with what is counted (<unit>, e.g.
   "reads"), with a legend giving the share of each part (the empty parts
   are left out; nothing is drawn if all are) */
void markup_bar(html_writer_t *w, const char *title, const char *unit,
    const markup_part_t *parts, size_t n);

/* A legend entry: a swatch of the class <cls> and its <label> */
void markup_swatch(html_writer_t *w, const char *cls, const char *label);

#endif
