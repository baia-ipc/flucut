/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include "html-page.h"
#include "html-report-markup.h"
#include "numutils.h"

void markup_open(html_writer_t *w, const char *tag, const char *cls) {
  html_open(w, tag);
  html_attr(w, "class", "%s", cls);
}

void markup_textf(html_writer_t *w, const char *tag, const char *cls,
    const char *fmt, ...) {
  markup_open(w, tag, cls);
  va_list ap;
  va_start(ap, fmt);
  html_vtextf(w, fmt, ap);
  va_end(ap);
  html_close(w);
}

void markup_note(html_writer_t *w, const char *fmt, ...) {
  markup_open(w, "p", "note");
  va_list ap;
  va_start(ap, fmt);
  html_vtextf(w, fmt, ap);
  va_end(ap);
  html_close(w);
}

void markup_open_row(html_writer_t *w) {
  html_open(w, "tr");
  html_inline(w);
}

void markup_header_row(html_writer_t *w, ...) {
  markup_open_row(w);
  va_list ap;
  va_start(ap, w);
  for (const char *h = va_arg(ap, const char *); h != NULL;
       h = va_arg(ap, const char *)) {
    html_open(w, "th");
    html_raw(w, h);
    html_close(w);
  }
  va_end(ap);
  html_close(w);
}

void markup_td_none(html_writer_t *w) {
  html_open(w, "td");
  html_raw(w, "&ndash;");
  html_close(w);
}

void markup_tile(html_writer_t *w, const char *label, const char *sub,
    const char *fmt, ...) {
  markup_open(w, "div", "tile");
  markup_open(w, "div", "v");
  va_list ap;
  va_start(ap, fmt);
  html_vtextf(w, fmt, ap);
  va_end(ap);
  html_close(w);
  markup_textf(w, "div", "l", "%s", label);
  if (sub != NULL)
    markup_textf(w, "div", "u", "%s", sub);
  html_close(w);
}

void markup_badge(html_writer_t *w, markup_state_t state, const char *fmt,
    ...) {
  /* (a check mark, an exclamation mark, a cross, a dot) */
  static const char *const cls[] = {"st ok", "st warn", "st bad", "st na"};
  static const char *const symbol[] = {"\xe2\x9c\x93", "!", "\xe2\x9c\x95",
                                       "\xc2\xb7"};
  char text[128];
  va_list ap;
  va_start(ap, fmt);
  vsnprintf(text, sizeof(text), fmt, ap);
  va_end(ap);
  markup_open(w, "span", cls[state]);
  html_textf(w, "%s%s%s", symbol[state], text[0] != '\0' ? " " : "", text);
  html_close(w);
}

void markup_details_start(html_writer_t *w, const char *summary) {
  html_open(w, "details");
  html_tag_text(w, "summary", summary);
}

void markup_details_end(html_writer_t *w) {
  html_close(w);
}

void markup_swatch(html_writer_t *w, const char *cls, const char *label) {
  html_open(w, "span");
  html_open(w, "i");
  html_attr(w, "class", "%s", cls);
  html_close(w);
  html_text(w, label);
  html_close(w);
}

void markup_bar(html_writer_t *w, const char *title, const char *unit,
    const markup_part_t *parts, size_t n) {
  assert(parts != NULL);
  size_t total = 0;
  for (size_t i = 0; i < n; i++)
    total += parts[i].n;
  if (total == 0)
    return;
  markup_open(w, "div", "bar");
  markup_open(w, "div", "bl");
  html_textf(w, "%s ", title);
  markup_textf(w, "span", "note", "of %s %s", html_num(total).s, unit);
  html_close(w);
  markup_open(w, "div", "sbar");
  for (size_t i = 0; i < n; i++) {
    if (parts[i].n == 0)
      continue;
    html_open(w, "i");
    html_attr(w, "class", "%s", parts[i].cls);
    html_attr(w, "style", "width:%.2f%%", html_percent(parts[i].n, total));
    html_attr(w, "title", "%s: %s (%s)", parts[i].label,
              html_num(parts[i].n).s, percent_text(parts[i].n, total).s);
    html_close(w);
  }
  html_close(w);
  markup_open(w, "div", "lg");
  for (size_t i = 0; i < n; i++) {
    if (parts[i].n == 0)
      continue;
    html_open(w, "span");
    html_attr(w, "title", "%s", html_num(parts[i].n).s);
    html_open(w, "i");
    html_attr(w, "class", "%s", parts[i].cls);
    html_close(w);
    html_textf(w, "%s ", parts[i].label);
    html_tag_text(w, "b", percent_text(parts[i].n, total).s);
    html_close(w);
  }
  html_close(w);
  html_close(w);
}
