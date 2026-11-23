/* (c) 2026 Giorgio Gonnella */

/*
   The elements shared by the plots of the HTML report (inline SVG): the
   lines of the axes and marks, and the labels.
*/

#ifndef HTML_REPORT_SVG_H
#define HTML_REPORT_SVG_H

#include "html-writer.h"

/* A line from (x1, y1) to (x2, y2), of the class <cls> (axis, grid,
   mark, site) */
void svg_line(html_writer_t *w, const char *cls, double x1, double y1,
    double x2, double y2);

/* A label at (x, y), anchored as <anchor> (start, middle, end), with the
   text formatted as by printf */
void svg_label(html_writer_t *w, double x, double y, const char *anchor,
    const char *fmt, ...) __attribute__((format(printf, 5, 6)));

#endif
