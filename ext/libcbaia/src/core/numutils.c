/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "numutils.h"
#include "test.h"

size_t hist_median(const size_t *hist, size_t n_bins, size_t bin_width,
    size_t total) {
  assert(hist != NULL || n_bins == 0);
  if (n_bins == 0)
    return 0;
  size_t cum = 0;
  for (size_t b = 0; b < n_bins; b++) {
    cum += hist[b];
    if (2 * cum >= total)
      return b * bin_width + bin_width / 2;
  }
  /* total exceeds the counts: the median is beyond the counted values,
     i.e. at least in the last bin (which usually counts also the larger
     values) */
  return (n_bins - 1) * bin_width + bin_width / 2;
}

double nice_number(double v) {
  if (v <= 1)
    return 1;
  double p = pow(10, floor(log10(v))), m = v / p;
  return (m <= 1 ? 1 : m <= 2 ? 2 : m <= 5 ? 5 : 10) * p;
}

/* the text of <percent>, of a part which is none, all or neither */
static percent_text_t percent_format(double percent, bool none, bool all) {
  percent_text_t out;
  if (!none && percent < 0.1)
    snprintf(out.s, sizeof(out.s), "<0.1%%");
  else if (!all && percent > 99.9 && percent <= 100.0)
    snprintf(out.s, sizeof(out.s), ">99.9%%");
  else
    snprintf(out.s, sizeof(out.s), "%.1f%%", percent);
  return out;
}

percent_text_t percent_text(size_t part, size_t total) {
  if (total == 0)
    return percent_format(0.0, true, false);
  return percent_format(100.0 * (double)part / (double)total, part == 0,
                        part == total);
}

percent_text_t fraction_percent_text(double fraction) {
  return percent_format(100.0 * fraction, fraction == 0.0, fraction == 1.0);
}

bool parse_double(const char *s, double *value) {
  char *end;
  errno = 0;
  *value = strtod(s, &end);
  return *s != '\0' && *end == '\0' && errno == 0 && isfinite(*value);
}

bool parse_percent(const char *s, double *value) {
  char buf[32];
  size_t len = strlen(s);
  if (len < 2 || len > sizeof(buf) || s[len - 1] != '%')
    return false;
  memcpy(buf, s, len - 1);
  buf[len - 1] = '\0';
  double percent;
  if (!parse_double(buf, &percent) || percent < 0 || percent >= 100)
    return false;
  *value = percent / 100;
  return true;
}

void numutils_test(void) {
  test_eq_str(percent_text(0, 2000).s, "0.0%");
  test_eq_str(percent_text(1, 2000).s, "<0.1%");
  test_eq_str(percent_text(1, 1000).s, "0.1%");
  test_eq_str(percent_text(423, 1000).s, "42.3%");
  test_eq_str(percent_text(999, 1000).s, "99.9%");
  test_eq_str(percent_text(1999, 2000).s, ">99.9%");
  test_eq_str(percent_text(2000, 2000).s, "100.0%");
  test_eq_str(percent_text(999999999, 1000000000).s, ">99.9%");
  test_eq_str(percent_text(3, 2).s, "150.0%");
  test_eq_str(percent_text(0, 0).s, "0.0%");
  test_eq_str(fraction_percent_text(0.0).s, "0.0%");
  test_eq_str(fraction_percent_text(0.0004).s, "<0.1%");
  test_eq_str(fraction_percent_text(0.031).s, "3.1%");
  test_eq_str(fraction_percent_text(0.9996).s, ">99.9%");
  test_eq_str(fraction_percent_text(1.0).s, "100.0%");
  test_eq(floordiv(7, 2), 3L, long, "%ld");
  test_eq(floordiv(-7, 2), -4L, long, "%ld");
  test_eq(floordiv(7, -2), -4L, long, "%ld");
  test_eq(floordiv(-6, 2), -3L, long, "%ld");
  test_eq(floordiv(-1, 2), -1L, long, "%ld");
  test_eq(floordiv(0, 5), 0L, long, "%ld");
  size_t hist[] = {0, 3, 1, 5, 1};
  test_eq_size(hist_median(hist, 5, 10, 10), 35);
  test_eq_size(hist_median(hist, 5, 10, 4), 15);
  test_eq_size(hist_median(hist, 5, 10, 0), 5);
  size_t empty[] = {0, 0};
  test_eq_size(hist_median(empty, 2, 10, 0), 5);
  /* total larger than the counts (sum 10): the last bin */
  test_eq_size(hist_median(hist, 5, 10, 30), 45);
  test_eq_size(hist_median(empty, 2, 10, 1), 15);
  test_eq_size(hist_median(NULL, 0, 10, 1), 0);
  test_assert(nice_number(0.3) == 1);
  test_assert(nice_number(1) == 1);
  test_assert(nice_number(1.5) == 2);
  test_assert(nice_number(3) == 5);
  test_assert(nice_number(7) == 10);
  test_assert(nice_number(130) == 200);
  test_assert(nice_number(4000) == 5000);
}
