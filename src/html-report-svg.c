/* (c) 2026 Giorgio Gonnella */

#include <stdarg.h>
#include "html-report-svg.h"

void svg_line(html_writer_t *w, const char *cls, double x1,
    double y1, double x2, double y2) {
  html_elem(w, "line");
  html_attr(w, "class", "%s", cls);
  html_attr(w, "x1", "%.1f", x1);
  html_attr(w, "y1", "%.1f", y1);
  html_attr(w, "x2", "%.1f", x2);
  html_attr(w, "y2", "%.1f", y2);
}

void svg_label(html_writer_t *w, double x, double y,
    const char *anchor, const char *fmt, ...) {
  html_open(w, "text");
  html_attr(w, "class", "muted");
  html_attr(w, "x", "%.1f", x);
  html_attr(w, "y", "%.1f", y);
  html_attr(w, "text-anchor", "%s", anchor);
  html_inline(w);
  va_list ap;
  va_start(ap, fmt);
  html_vtextf(w, fmt, ap);
  va_end(ap);
  html_close(w);
}
