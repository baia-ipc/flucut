/* (c) 2026 Giorgio Gonnella */

/*
   HTML templates: static markup (e.g. the skeleton of a page) written as
   a normal, indented HTML file, embedded in the program at compile time
   (dev/embed.sh), with named placeholders filled when it is written:

     {{name}}      the value, escaped
     {{{name}}}    the value, not escaped (trusted markup, e.g. a style)
     {{#name}}     the following, up to {{/name}}, only if the value is
     {{/name}}     not NULL (a section; sections can be nested)

   Names: lower case letters, digits and _. A line containing only a
   section tag (and spaces) is not written. The template is written by a
   html_writer_t, indented as the current element; the lines of a value
   are indented as the line of its placeholder, and a newline at the end
   of a value is dropped, so that e.g. a style can be written as

     <style>
       {{{css}}}
     </style>

   A placeholder without value, a value not used by the template, a NULL
   value outside of a section of its name and a syntax error are
   programming errors: the program is aborted with a message, also when
   compiled with NDEBUG (the templates are constant: any test writing a
   template finds them).
*/

#ifndef HTML_TEMPLATE_H
#define HTML_TEMPLATE_H

#include <stddef.h>
#include "html-writer.h"

/* maximal number of values of a template and of nested sections */
#define HTML_TEMPLATE_MAX_VALUES 32
#define HTML_TEMPLATE_MAX_NESTING 8

/*
   Write the template <tmpl> filled with the values, given as name and
   value pairs of strings, e.g.:
     html_template(w, page_tmpl, "title", title, "css", css);
   (a macro: the pairs are counted at compile time)
*/
#define html_template(W, TMPL, ...) \
  html_template_values((W), (TMPL), \
      (const char *const[]){__VA_ARGS__}, \
      sizeof((const char *const[]){__VA_ARGS__}) / sizeof(const char *))

/* The same, with the <n> names and values in <kv> (name, value, name,
   value...) */
void html_template_values(html_writer_t *w, const char *tmpl,
    const char *const *kv, size_t n);

/* Unit tests */
void html_template_test(void);

#endif
