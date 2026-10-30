/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "html-writer.h"
#include "strutils.h"
#include "test.h"

/* content of an open element: none yet, elements on their own lines,
   or text (and then everything on the same line) */
enum { CONTENT_NONE, CONTENT_LINES, CONTENT_INLINE };

#define INDENT 2

void html_writer_init(html_writer_t *w, FILE *f) {
  assert(w != NULL);
  assert(f != NULL);
  memset(w, 0, sizeof(*w));
  w->f = f;
}

/* within an element written on one line */
static bool in_line(const html_writer_t *w) {
  return w->inline_from > 0;
}

static void indent(html_writer_t *w, size_t depth) {
  fprintf(w->f, "%*s", (int)(depth * INDENT), "");
}

/* complete the open start tag, if any */
static void end_start_tag(html_writer_t *w) {
  if (w->start == NULL)
    return;
  if (w->start_void) {
    /* (SVG is XML: an element without content is closed by "/>") */
    fputs(w->svg > 0 ? "/>" : ">", w->f);
    if (!in_line(w))
      fputc('\n', w->f);
  } else
    fputc('>', w->f);
  w->start = NULL;
}

size_t html_begin_lines(html_writer_t *w) {
  assert(w != NULL);
  end_start_tag(w);
  assert(!in_line(w));
  if (w->depth > 0) {
    unsigned char *c = &w->content[w->depth - 1];
    if (*c == CONTENT_NONE)
      fputc('\n', w->f);  /* the line of its start tag ends */
    *c = CONTENT_LINES;
  }
  return w->depth * INDENT;
}

/* start a line of content of the current element (outside of a line) */
static void begin_line(html_writer_t *w) {
  fprintf(w->f, "%*s", (int)html_begin_lines(w), "");
}

static void start_tag(html_writer_t *w, const char *tag, bool is_void) {
  assert(w != NULL);
  assert(tag != NULL);
  end_start_tag(w);
  if (!in_line(w))
    begin_line(w);
  fprintf(w->f, "<%s", tag);
  w->start = tag;
  w->start_void = is_void;
}

void html_open(html_writer_t *w, const char *tag) {
  assert(w != NULL);
  assert(tag != NULL);
  /* the open elements are stored in arrays of fixed size */
  if (w->depth >= HTML_WRITER_MAX_DEPTH)
    fatal_printf("html_open: <%s> nested too deeply (more than %d open "
                 "elements)\n", tag, HTML_WRITER_MAX_DEPTH);
  start_tag(w, tag, false);
  w->tags[w->depth] = tag;
  w->content[w->depth] = CONTENT_NONE;
  w->depth++;
  if (strcmp(tag, "svg") == 0)
    w->svg++;
}

void html_elem(html_writer_t *w, const char *tag) {
  start_tag(w, tag, true);
}

/* fputs_html_escaped() of the first <len> characters */
static void write_escaped(FILE *f, const char *s, size_t len) {
  for (size_t i = 0; i < len; i++)
    switch (s[i]) {
      case '&': fputs("&amp;", f); break;
      case '<': fputs("&lt;", f); break;
      case '>': fputs("&gt;", f); break;
      case '"': fputs("&quot;", f); break;
      default: fputc(s[i], f);
    }
}

/* write the formatted string, escaped (in a local buffer if it fits) */
static void write_escaped_vf(FILE *f, const char *fmt, va_list ap) {
  char local[256], *s = local;
  va_list ap2;
  va_copy(ap2, ap);
  int len = vsnprintf(local, sizeof(local), fmt, ap);
  assert(len >= 0);
  if ((size_t)len >= sizeof(local)) {
    XMALLOC(s, (size_t)len + 1);
    vsnprintf(s, (size_t)len + 1, fmt, ap2);
  }
  va_end(ap2);
  write_escaped(f, s, (size_t)len);
  if (s != local)
    free(s);
}

void html_attr(html_writer_t *w, const char *name, const char *fmt, ...) {
  assert(w != NULL);
  assert(name != NULL);
  assert(fmt != NULL);
  assert(w->start != NULL);  /* only directly after the start of a tag */
  fprintf(w->f, " %s=\"", name);
  va_list ap;
  va_start(ap, fmt);
  write_escaped_vf(w->f, fmt, ap);
  va_end(ap);
  fputc('"', w->f);
}

void html_inline(html_writer_t *w) {
  assert(w != NULL);
  assert(w->depth > 0);
  assert(w->content[w->depth - 1] == CONTENT_NONE);
  w->content[w->depth - 1] = CONTENT_INLINE;
  if (!in_line(w))
    w->inline_from = w->depth;
}

/* prepare the writing of text in the current element: on the line of its
   start tag, or on a line of its own after elements on their own lines;
   true in the second case (the line is then ended by the caller) */
static bool begin_text(html_writer_t *w) {
  end_start_tag(w);
  if (in_line(w))
    return false;
  if (w->depth > 0 && w->content[w->depth - 1] == CONTENT_NONE) {
    html_inline(w);
    return false;
  }
  begin_line(w);
  return true;
}

void html_text(html_writer_t *w, const char *text) {
  assert(w != NULL);
  assert(text != NULL);
  bool own_line = begin_text(w);
  write_escaped(w->f, text, strlen(text));
  if (own_line)
    fputc('\n', w->f);
}

void html_vtextf(html_writer_t *w, const char *fmt, va_list ap) {
  assert(w != NULL);
  assert(fmt != NULL);
  bool own_line = begin_text(w);
  write_escaped_vf(w->f, fmt, ap);
  if (own_line)
    fputc('\n', w->f);
}

void html_textf(html_writer_t *w, const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  html_vtextf(w, fmt, ap);
  va_end(ap);
}

void html_raw(html_writer_t *w, const char *markup) {
  assert(w != NULL);
  assert(markup != NULL);
  bool own_line = begin_text(w);
  fputs(markup, w->f);
  if (own_line)
    fputc('\n', w->f);
}

void html_close(html_writer_t *w) {
  assert(w != NULL);
  assert(w->depth > 0);
  end_start_tag(w);
  w->depth--;
  const char *tag = w->tags[w->depth];
  if (w->content[w->depth] == CONTENT_LINES)
    indent(w, w->depth);
  fprintf(w->f, "</%s>", tag);
  if (w->inline_from > w->depth)
    w->inline_from = 0;
  if (!in_line(w))
    fputc('\n', w->f);
  if (strcmp(tag, "svg") == 0)
    w->svg--;
}

void html_tag_text(html_writer_t *w, const char *tag, const char *text) {
  html_open(w, tag);
  html_text(w, text);
  html_close(w);
}

void html_tag_textf(html_writer_t *w, const char *tag, const char *fmt,
    ...) {
  html_open(w, tag);
  va_list ap;
  va_start(ap, fmt);
  html_vtextf(w, fmt, ap);
  va_end(ap);
  html_close(w);
}

void html_writer_end(html_writer_t *w) {
  assert(w != NULL);
  end_start_tag(w);
  assert(w->depth == 0);
}

/* --- unit tests --- */

/* the output of the writer, in <out> (to be freed) */
#define WRITER_TEST_START \
  char *out = NULL; \
  size_t size = 0; \
  FILE *f = open_memstream(&out, &size); \
  test_assert(f != NULL); \
  html_writer_t w; \
  html_writer_init(&w, f)

#define WRITER_TEST_END(EXP) do { \
  html_writer_end(&w); \
  fclose(f); \
  test_eq_str(out, EXP); \
  free(out); \
} while (0)

void html_writer_test(void) {
  /* nesting, attributes, escaping */
  {
    WRITER_TEST_START;
    html_open(&w, "div");
    html_attr(&w, "class", "%s", "a\"b");
    html_open(&w, "p");
    html_text(&w, "x < y & z");
    html_close(&w);
    html_open(&w, "table");
    html_close(&w);
    html_close(&w);
    WRITER_TEST_END("<div class=\"a&quot;b\">\n"
                    "  <p>x &lt; y &amp; z</p>\n"
                    "  <table></table>\n"
                    "</div>\n");
  }
  /* SVG: elements without content, formatted attributes and text */
  {
    WRITER_TEST_START;
    html_open(&w, "svg");
    html_attr(&w, "viewBox", "0 0 %d %d", 760, 50);
    html_elem(&w, "rect");
    html_attr(&w, "x", "%.1f", 2.25);
    html_attr(&w, "y", "%d", 3);
    html_open(&w, "rect");
    html_attr(&w, "class", "s");
    html_tag_textf(&w, "title", "%s: %d", "a", 1);
    html_close(&w);
    html_open(&w, "text");
    html_textf(&w, "%d%% <", 50);
    html_close(&w);
    html_close(&w);
    WRITER_TEST_END("<svg viewBox=\"0 0 760 50\">\n"
                    "  <rect x=\"2.2\" y=\"3\"/>\n"
                    "  <rect class=\"s\">\n"
                    "    <title>a: 1</title>\n"
                    "  </rect>\n"
                    "  <text>50% &lt;</text>\n"
                    "</svg>\n");
  }
  /* HTML void elements; mixed content and inline elements on one line */
  {
    WRITER_TEST_START;
    html_elem(&w, "meta");
    html_attr(&w, "charset", "utf-8");
    html_open(&w, "p");
    html_text(&w, "a ");
    html_tag_text(&w, "b", "b");
    html_raw(&w, " &middot; ");
    html_open(&w, "span");
    html_attr(&w, "class", "note");
    html_elem(&w, "img");
    html_attr(&w, "src", "x.png");
    html_close(&w);
    html_close(&w);
    html_open(&w, "tr");
    html_inline(&w);
    html_tag_text(&w, "td", "1");
    html_tag_text(&w, "td", "2");
    html_close(&w);
    WRITER_TEST_END("<meta charset=\"utf-8\">\n"
                    "<p>a <b>b</b> &middot; <span class=\"note\">"
                    "<img src=\"x.png\"></span></p>\n"
                    "<tr><td>1</td><td>2</td></tr>\n");
  }
  /* text after elements on their own lines; long attribute values */
  {
    char long_value[1000];
    memset(long_value, 'x', sizeof(long_value) - 1);
    long_value[sizeof(long_value) - 1] = '\0';
    WRITER_TEST_START;
    html_open(&w, "div");
    html_tag_text(&w, "h3", "t");
    html_text(&w, "after");
    html_elem(&w, "img");
    html_attr(&w, "src", "%s<", long_value);
    html_close(&w);
    html_writer_end(&w);
    fclose(f);
    const char *exp_start = "<div>\n  <h3>t</h3>\n  after\n  <img src=\"",
               *exp_end = "&lt;\">\n</div>\n";
    test_eq_size(strlen(out), strlen(exp_start) + sizeof(long_value) - 1 +
                 strlen(exp_end));
    test_assert(strncmp(out, exp_start, strlen(exp_start)) == 0);
    test_eq_str(out + strlen(out) - strlen(exp_end), exp_end);
    free(out);
  }
  /* the maximum number of open elements */
  {
    WRITER_TEST_START;
    for (size_t i = 0; i < HTML_WRITER_MAX_DEPTH; i++)
      html_open(&w, "div");
    test_eq_size(w.depth, (size_t)HTML_WRITER_MAX_DEPTH);
    for (size_t i = 0; i < HTML_WRITER_MAX_DEPTH; i++)
      html_close(&w);
    html_writer_end(&w);
    fclose(f);
    /* the innermost element (empty) on one line */
    char *exp = NULL;
    size_t exp_size = 0;
    FILE *e = open_memstream(&exp, &exp_size);
    test_assert(e != NULL);
    for (int i = 0; i < HTML_WRITER_MAX_DEPTH - 1; i++)
      fprintf(e, "%*s<div>\n", INDENT * i, "");
    fprintf(e, "%*s<div></div>\n", INDENT * (HTML_WRITER_MAX_DEPTH - 1), "");
    for (int i = HTML_WRITER_MAX_DEPTH - 2; i >= 0; i--)
      fprintf(e, "%*s</div>\n", INDENT * i, "");
    fclose(e);
    test_eq_str(out, exp);
    free(exp);
    free(out);
  }
}
