/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "html-template.h"
#include "test.h"

/* --- parsing --- */

typedef enum {
  TOK_END, TOK_TEXT, TOK_VALUE, TOK_RAW, TOK_OPEN, TOK_CLOSE, TOK_ERROR
} token_type_t;

typedef struct {
  token_type_t type;
  const char *s;  /* the text, or the name of the placeholder */
  size_t len;
} token_t;

typedef struct {
  const char *pos;
  bool line_start;
} scanner_t;

static bool is_name_char(char c) {
  return (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_';
}

/* the tag starting at <p> ("{{"), in <t>, and the position after it in
   <end>; false if it is not valid */
static bool scan_tag(const char *p, token_t *t, const char **end) {
  bool raw = strncmp(p, "{{{", 3) == 0;
  const char *close = raw ? "}}}" : "}}";
  p += raw ? 3 : 2;
  t->type = raw ? TOK_RAW : TOK_VALUE;
  if (!raw && (*p == '#' || *p == '/')) {
    t->type = *p == '#' ? TOK_OPEN : TOK_CLOSE;
    p++;
  }
  t->s = p;
  while (is_name_char(*p))
    p++;
  t->len = (size_t)(p - t->s);
  if (t->len == 0 || strncmp(p, close, strlen(close)) != 0)
    return false;
  *end = p + strlen(close);
  return true;
}

/* the next token: a placeholder, a section tag, or text up to the next
   tag or the end of the line (included) */
static void next_token(scanner_t *sc, token_t *t) {
  const char *p = sc->pos, *end;
  if (*p == '\0') {
    t->type = TOK_END;
    return;
  }
  /* a line containing only a section tag is skipped */
  if (sc->line_start) {
    const char *q = p + strspn(p, " ");
    if (strncmp(q, "{{", 2) == 0 && (q[2] == '#' || q[2] == '/') &&
        scan_tag(q, t, &end)) {
      end += strspn(end, " ");
      if (*end == '\n' || *end == '\0') {
        sc->pos = *end == '\n' ? end + 1 : end;
        return;
      }
    }
  }
  if (strncmp(p, "{{", 2) == 0) {
    if (!scan_tag(p, t, &end)) {
      t->type = TOK_ERROR;
      t->s = p;
      t->len = strcspn(p, "\n");
      return;
    }
    sc->pos = end;
    sc->line_start = false;
    return;
  }
  end = p;
  while (*end != '\0' && *end != '\n' && strncmp(end, "{{", 2) != 0)
    end++;
  if (*end == '\n')
    end++;
  t->type = TOK_TEXT;
  t->s = p;
  t->len = (size_t)(end - p);
  sc->pos = end;
  sc->line_start = end[-1] == '\n';
}

/* --- check --- */

typedef struct {
  const char *msg;
  const char *at;  /* where (in the template, or the name of a value) */
  size_t len;
} template_error_t;

/* index of the value named as the token, or -1 */
static long find_value(const char *const *kv, size_t n, const token_t *t) {
  for (size_t i = 0; i < n; i += 2)
    if (strlen(kv[i]) == t->len && strncmp(kv[i], t->s, t->len) == 0)
      return (long)(i / 2);
  return -1;
}

static bool set_error(template_error_t *err, const char *msg,
    const char *at, size_t len) {
  err->msg = msg;
  err->at = at;
  err->len = len;
  return false;
}

/* check the template and the values (see html-template.h); false (and
   the error in <err>) if not valid */
static bool template_check(const char *tmpl, const char *const *kv,
    size_t n, template_error_t *err) {
  if (n % 2 != 0 || n / 2 > HTML_TEMPLATE_MAX_VALUES)
    return set_error(err, "odd or too large number of names and values",
                     "", 0);
  for (size_t i = 0; i < n; i += 2) {
    if (kv[i] == NULL)
      return set_error(err, "NULL name", "", 0);
    token_t t = {TOK_VALUE, kv[i], strlen(kv[i])};
    if (find_value(kv, i, &t) >= 0)
      return set_error(err, "name given twice", kv[i], t.len);
  }
  bool used[HTML_TEMPLATE_MAX_VALUES] = {false};
  long open[HTML_TEMPLATE_MAX_NESTING];
  size_t depth = 0, skip = 0;
  scanner_t sc = {tmpl, true};
  token_t t;
  for (next_token(&sc, &t); t.type != TOK_END; next_token(&sc, &t)) {
    if (t.type == TOK_TEXT)
      continue;
    if (t.type == TOK_ERROR)
      return set_error(err, "invalid tag", t.s, t.len);
    long v = find_value(kv, n, &t);
    if (v < 0)
      return set_error(err, "no value for the placeholder", t.s, t.len);
    used[v] = true;
    bool null = kv[2 * v + 1] == NULL;
    switch (t.type) {
      case TOK_VALUE:
      case TOK_RAW:
        if (null && skip == 0)
          return set_error(err, "NULL value outside of its section", t.s,
                           t.len);
        break;
      case TOK_OPEN:
        if (depth == HTML_TEMPLATE_MAX_NESTING)
          return set_error(err, "too many nested sections", t.s, t.len);
        open[depth++] = v;
        skip += null;
        break;
      default:  /* TOK_CLOSE */
        if (depth == 0 || open[depth - 1] != v)
          return set_error(err, "end of a section not open", t.s, t.len);
        depth--;
        skip -= null;
    }
  }
  if (depth > 0)
    return set_error(err, "section not closed", kv[2 * open[depth - 1]],
                     strlen(kv[2 * open[depth - 1]]));
  for (size_t i = 0; i < n / 2; i++)
    if (!used[i])
      return set_error(err, "value not used", kv[2 * i], strlen(kv[2 * i]));
  size_t len = strlen(tmpl);
  if (len == 0 || tmpl[len - 1] != '\n')
    return set_error(err, "the template does not end with a newline",
                     "", 0);
  return true;
}

/* --- filling --- */

/* output: the indentation of the writer is written before each line, and
   for the lines of a value also the indentation of the line of its
   placeholder (only before a character: no trailing spaces) */
typedef struct {
  FILE *f;
  size_t indent;   /* of the writer */
  size_t pending;  /* spaces to write before the next character */
  size_t lead;     /* spaces at the start of the current template line */
  bool in_lead;
} output_t;

static void put(output_t *o, char c, bool from_template, bool escape) {
  if (c == '\n') {
    fputc('\n', o->f);
    o->pending = o->indent + (from_template ? 0 : o->lead);
    if (from_template) {
      o->lead = 0;
      o->in_lead = true;
    }
    return;
  }
  if (o->pending > 0) {
    fprintf(o->f, "%*s", (int)o->pending, "");
    o->pending = 0;
  }
  if (from_template && o->in_lead && c == ' ')
    o->lead++;
  else
    o->in_lead = false;
  if (!escape) {
    fputc(c, o->f);
    return;
  }
  switch (c) {
    case '&': fputs("&amp;", o->f); break;
    case '<': fputs("&lt;", o->f); break;
    case '>': fputs("&gt;", o->f); break;
    case '"': fputs("&quot;", o->f); break;
    default: fputc(c, o->f);
  }
}

void html_template_values(html_writer_t *w, const char *tmpl,
    const char *const *kv, size_t n) {
  assert(w != NULL);
  assert(tmpl != NULL);
  assert(kv != NULL || n == 0);
  template_error_t err;
  if (!template_check(tmpl, kv, n, &err)) {
    fprintf(stderr, "html_template: %s: \"%.*s\"\n", err.msg,
            (int)err.len, err.at);
    abort();
  }
  output_t o = {w->f, html_begin_lines(w), 0, 0, true};
  o.pending = o.indent;
  bool null[HTML_TEMPLATE_MAX_NESTING];
  size_t depth = 0, skip = 0;
  scanner_t sc = {tmpl, true};
  token_t t;
  for (next_token(&sc, &t); t.type != TOK_END; next_token(&sc, &t)) {
    const char *value = NULL;
    if (t.type != TOK_TEXT)
      value = kv[2 * find_value(kv, n, &t) + 1];
    switch (t.type) {
      case TOK_TEXT:
        if (skip == 0)
          for (size_t i = 0; i < t.len; i++)
            put(&o, t.s[i], true, false);
        break;
      case TOK_VALUE:
      case TOK_RAW:
        if (skip == 0) {
          size_t len = strlen(value);
          if (len > 0 && value[len - 1] == '\n')
            len--;
          for (size_t i = 0; i < len; i++)
            put(&o, value[i], false, t.type == TOK_VALUE);
        }
        break;
      case TOK_OPEN:
        null[depth++] = value == NULL;
        skip += value == NULL;
        break;
      default:  /* TOK_CLOSE (the template was checked) */
        skip -= null[--depth];
    }
  }
}

/* --- unit tests --- */

/* template_check() of <tmpl> with the values fails with <msg> */
#define TEMPLATE_CHECK_FAILS(TMPL, MSG, ...) do { \
  const char *const kv[] = {__VA_ARGS__}; \
  template_error_t err; \
  test_assert(!template_check(TMPL, kv, sizeof(kv) / sizeof(kv[0]), \
                              &err)); \
  test_eq_str(err.msg, MSG); \
} while (0)

void html_template_test(void) {
  char *out = NULL;
  size_t size = 0;
  FILE *f = open_memstream(&out, &size);
  test_assert(f != NULL);
  html_writer_t w;
  html_writer_init(&w, f);
  html_open(&w, "div");
  html_template(&w, "<p>{{a}} {{{b}}}</p>\n"
                    "<div>\n"
                    "  {{{c}}}\n"
                    "</div>\n"
                    "{{#s}}\n"
                    "<p>{{s}}</p>\n"
                    "{{/s}}\n"
                    "  {{#n}}\n"
                    "<p>{{n}}</p>\n"
                    "{{/n}}\n"
                    "<p>{{#s}}{{s}}{{/s}}{{#n}}x{{/n}}.</p>\n",
                "a", "x<y", "b", "<b>&amp;</b>",
                "c", "<i>1</i>\n\n<i>2</i>\n", "s", "\"s\"", "n", NULL);
  html_close(&w);
  html_writer_end(&w);
  fclose(f);
  test_eq_str(out, "<div>\n"
                   "  <p>x&lt;y <b>&amp;</b></p>\n"
                   "  <div>\n"
                   "    <i>1</i>\n"
                   "\n"
                   "    <i>2</i>\n"
                   "  </div>\n"
                   "  <p>&quot;s&quot;</p>\n"
                   "  <p>&quot;s&quot;.</p>\n"
                   "</div>\n");
  free(out);
  /* errors */
  template_error_t err;
  test_assert(template_check("<p>{{a}}</p>\n", (const char *[]){"a", "1"},
                             2, &err));
  TEMPLATE_CHECK_FAILS("<p>{{a}}</p>\n", "no value for the placeholder",
                       "b", "1");
  TEMPLATE_CHECK_FAILS("<p>{{a}}</p>\n", "value not used",
                       "a", "1", "b", "2");
  TEMPLATE_CHECK_FAILS("<p>{{a}}</p>\n", "NULL value outside of its section",
                       "a", NULL);
  TEMPLATE_CHECK_FAILS("{{#a}}\n<p>\n", "section not closed", "a", "1");
  TEMPLATE_CHECK_FAILS("{{#a}}{{#b}}{{/a}}{{/b}}\n",
                       "end of a section not open", "a", "1", "b", "2");
  TEMPLATE_CHECK_FAILS("<p>{{a}</p>\n", "invalid tag", "a", "1");
  TEMPLATE_CHECK_FAILS("<p>{{ a }}</p>\n", "invalid tag", "a", "1");
  TEMPLATE_CHECK_FAILS("<p>{{{a}}</p>\n", "invalid tag", "a", "1");
  TEMPLATE_CHECK_FAILS("<p>{{A}}</p>\n", "invalid tag", "A", "1");
  TEMPLATE_CHECK_FAILS("<p>{{a}}</p>\n", "name given twice",
                       "a", "1", "a", "2");
  TEMPLATE_CHECK_FAILS("<p>{{a}}</p>\n",
                       "odd or too large number of names and values",
                       "a", "1", "b");
  TEMPLATE_CHECK_FAILS("<p>{{a}}</p>",
                       "the template does not end with a newline", "a", "1");
}
