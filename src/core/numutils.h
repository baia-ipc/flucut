/* (c) 2026 Giorgio Gonnella */

/*
   Numeric utilities: integer division rounding towards minus infinity, the
   median of a histogram, the rounding of a value to a "nice" number
   (for the axes of plots), and the text of a percentage.
*/

#ifndef NUMUTILS_H
#define NUMUTILS_H

#include <stddef.h>
#include <stdbool.h>

/* a / b rounded towards minus infinity (as // in Python); b != 0 */
static inline long floordiv(long a, long b) {
  long q = a / b;
  return (a % b != 0 && ((a < 0) != (b < 0))) ? q - 1 : q;
}

/*
   Median of the values counted in a histogram of <n_bins> bins of width
   <bin_width> (bin b: values from b * bin_width to (b + 1) * bin_width - 1),
   with <total> values in all: the middle of the bin containing the median
   (the middle of the first bin if total is 0); if total exceeds the sum of
   the counts and the median is not among the counted values, the middle
   of the last bin (usually counting also the larger values); 0 if n_bins
   is 0
*/
size_t hist_median(const size_t *hist, size_t n_bins, size_t bin_width,
    size_t total);

/*
   The smallest "nice" number not below <v>, i.e. 1, 2 or 5 times a power
   of 10 (e.g. for the maximum or the step of the axis of a plot):
   0.3 -> 1, 1.5 -> 2, 3 -> 5, 7 -> 10, 130 -> 200; values up to 1 give 1
   (the result is never below 1)
*/
double nice_number(double v);

/* a percentage formatted by percent_text() or fraction_percent_text() */
typedef struct {
  char s[16];
} percent_text_t;

/*
   The percentage <part> of <total> with one decimal, e.g. "42.3%": the
   string is percent_text(a, b).s, valid up to the end of the statement (a
   value, not a static buffer). "0.0%" and "100.0%" are written only for
   none and for all: a part which is not zero, but below 0.1%, gives
   "<0.1%", and one above 99.9%, but not all, gives ">99.9%". If total is
   0, "0.0%"
*/
percent_text_t percent_text(size_t part, size_t total);

/* The same, for a <fraction> of 1 (0.0: none, 1.0: all) */
percent_text_t fraction_percent_text(double fraction);

/* Parse a finite decimal number written in the whole string. */
bool parse_double(const char *s, double *value);

/* Parse a percentage in [0, 100), written with a trailing '%', as a
   fraction in [0, 1). */
bool parse_percent(const char *s, double *value);

/* Unit tests */
void numutils_test(void);

#endif
