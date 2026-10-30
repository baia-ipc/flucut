/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <math.h>
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

void numutils_test(void) {
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
