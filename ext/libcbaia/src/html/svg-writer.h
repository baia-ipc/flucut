/* (c) 2026 Giorgio Gonnella */

/* Small SVG primitives written through the HTML writer. */

#ifndef SVG_WRITER_H
#define SVG_WRITER_H

#include "html-writer.h"

/* A line from (x1, y1) to (x2, y2), of the class <cls>. */
void svg_line(html_writer_t *w, const char *cls, double x1, double y1,
    double x2, double y2);

/* A label at (x, y), anchored as <anchor> (start, middle, end), with the
   text formatted as by printf. */
void svg_label(html_writer_t *w, double x, double y, const char *anchor,
    const char *fmt, ...) __attribute__((format(printf, 5, 6)));

#endif
