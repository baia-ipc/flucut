/* (c) 2026 Giorgio Gonnella */

/*
   Writing HTML (and inline SVG) element by element: the text and the
   attribute values are escaped in one place, and the output is indented,
   one element per line, so that the nesting is visible both in the code
   and in the page.

     html_writer_t w;                     output:
     html_writer_init(&w, f);
     html_open(&w, "svg");                <svg viewBox="0 0 760 50">
     html_attr(&w, "viewBox", "0 0 %d %d", 760, 50);
     html_elem(&w, "rect");                 <rect class="s" width="20"/>
     html_attr(&w, "class", "%s", "s");
     html_attr(&w, "width", "%d", 20);
     html_open(&w, "text");                 <text x="4">a &lt; b</text>
     html_attr(&w, "x", "%d", 4);
     html_text(&w, "a < b");
     html_close(&w);
     html_close(&w);                      </svg>
     html_writer_end(&w);

   Attributes: after html_open() or html_elem() the start tag stays open,
   and each attribute is added by a html_attr() call, whose printf format
   is checked by the compiler (as it would not be for name/value pairs in
   variable arguments); the start tag is completed by the next call of any
   other function.

   Layout: each element starts on a new line, indented by 2 spaces for
   each enclosing element, and its content follows on the next lines;
   except that an element containing text (e.g. p, td, h2, SVG text) is
   written on one line together with all its content, e.g. <p>a
   <b>b</b></p>, as line breaks between text and inline elements would
   change the rendering; html_inline() does the same for an element
   without text (e.g. a table row).

   The tag and attribute names must be string literals (or live as long
   as the element is open); at most HTML_WRITER_MAX_DEPTH elements can be
   open at the same time (opening one more exits the program with an error
   message).
*/

#ifndef HTML_WRITER_H
#define HTML_WRITER_H

#include <stdarg.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

#define HTML_WRITER_MAX_DEPTH 32

typedef struct {
  FILE *f;
  const char *tags[HTML_WRITER_MAX_DEPTH];  /* the open elements */
  unsigned char content[HTML_WRITER_MAX_DEPTH];  /* their content so far */
  size_t depth;
  size_t inline_from;  /* depth of the outermost element written on one
                          line (0: none) */
  size_t svg;          /* open svg elements */
  const char *start;   /* element whose start tag is open, or NULL */
  bool start_void;
} html_writer_t;

/* Initialize a writer to <f> */
void html_writer_init(html_writer_t *w, FILE *f);

/* Complete the output (all elements must be closed) */
void html_writer_end(html_writer_t *w);

/* Start an element (its attributes follow, then its content); end it
   with html_close() */
void html_open(html_writer_t *w, const char *tag);

/* An element without content nor end tag (its attributes follow): an
   HTML void element (e.g. img, written as <img ...>) or, in an svg
   element, any element (<rect .../>) */
void html_elem(html_writer_t *w, const char *tag);

/* An attribute of the element just started: its value is formatted as by
   printf() and escaped */
void html_attr(html_writer_t *w, const char *name, const char *fmt, ...)
  __attribute__((format(printf, 3, 4)));

/* Write the content of the current element on the same line as its tags
   (before any content) */
void html_inline(html_writer_t *w);

/* Text (escaped) in the current element */
void html_text(html_writer_t *w, const char *text);

/* Text formatted as by printf() (escaped) in the current element */
void html_textf(html_writer_t *w, const char *fmt, ...)
  __attribute__((format(printf, 2, 3)));

/* The same, with the arguments in a va_list */
void html_vtextf(html_writer_t *w, const char *fmt, va_list ap);

/* Trusted markup (not escaped), written as text */
void html_raw(html_writer_t *w, const char *markup);

/* End the innermost open element */
void html_close(html_writer_t *w);

/* An element without attributes, containing only <text> (escaped),
   e.g. <h2>Results</h2> */
void html_tag_text(html_writer_t *w, const char *tag, const char *text);

/* The same, with the text formatted as by printf() */
void html_tag_textf(html_writer_t *w, const char *tag, const char *fmt, ...)
  __attribute__((format(printf, 3, 4)));

/*
   Prepare the writing of lines of markup (e.g. a template, see
   html-template.h) as content of the current element: returns the
   indentation of the lines (number of spaces); the lines are then written
   to w->f, each one ending with a newline
*/
size_t html_begin_lines(html_writer_t *w);

/* Unit tests */
void html_writer_test(void);

#endif
