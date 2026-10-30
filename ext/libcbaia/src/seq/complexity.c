/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <math.h>
#include "complexity.h"
#include "test.h"

/* distinct molecules expected in n reads of a library of m molecules */
static double expected_distinct(double m, double n) {
  return m * -expm1(-n / m);
}

void complexity_estimate(size_t reads, double distinct_keys,
    double mean_key_length, double error_rate, complexity_t *out) {
  assert(out != NULL);
  assert(error_rate >= 0 && error_rate < 1);
  double n = (double)reads;
  out->error_free = pow(1.0 - error_rate, mean_key_length);
  double n_correct = out->error_free * n;
  /* (too few reads takes precedence) */
  out->high_error_rate = reads >= COMPLEXITY_MIN_READS &&
                         out->error_free < COMPLEXITY_MIN_ERROR_FREE;
  out->estimated = n_correct >= COMPLEXITY_MIN_READS && !out->high_error_rate;
  out->saturated_unknown = false;
  out->molecules_seen = out->duplication = out->molecules =
    out->molecules_2x = 0;
  if (!out->estimated)
    return;
  double d = distinct_keys - (1.0 - out->error_free) * n;
  if (d < 1)
    d = 1;
  if (d > n_correct)
    d = n_correct;
  /* no duplicates (within the error of the count): M is unknown */
  if (1.0 - d / n_correct < 0.005) {
    out->saturated_unknown = true;
    out->molecules_seen = n;
    return;
  }
  /* solve expected_distinct(m, n_correct) = d: increasing in m */
  double lo = d, hi = d;
  while (expected_distinct(hi, n_correct) < d)
    hi *= 2;
  for (int i = 0; i < 200 && hi - lo > 1e-6 * hi; i++) {
    double mid = (lo + hi) / 2;
    if (expected_distinct(mid, n_correct) < d)
      lo = mid;
    else
      hi = mid;
  }
  out->molecules = (lo + hi) / 2;
  /* the distinct molecules and duplicates in all reads */
  out->molecules_seen = expected_distinct(out->molecules, n);
  out->duplication = 1.0 - out->molecules_seen / n;
  out->molecules_2x = expected_distinct(out->molecules, 2 * n);
}

/* --- unit tests --- */

void complexity_test(void) {
  complexity_t c;
  /* too few reads */
  complexity_estimate(50, 50, 24, 0.0, &c);
  test_assert(!c.estimated);
  /* no errors: 10000 molecules, 20000 reads: 8647 distinct expected */
  complexity_estimate(20000, 8646.6, 24, 0.0, &c);
  test_assert(c.estimated && !c.saturated_unknown);
  test_assert(fabs(c.molecules - 10000) < 10);
  test_assert(fabs(c.duplication - (1 - 8646.6 / 20000)) < 1e-3);
  test_assert(fabs(c.molecules_2x - 9816.8) < 10);
  /* errors: p = 0.97^24 = 0.4814; 9628 correct reads of 10000 molecules
     give 6181 distinct; the 10372 reads with errors add 10372 keys */
  double p = pow(0.97, 24);
  double n_correct = p * 20000;
  double d = 10000 * -expm1(-n_correct / 10000);
  complexity_estimate(20000, d + (1 - p) * 20000, 24, 0.03, &c);
  test_assert(fabs(c.error_free - p) < 1e-9);
  test_assert(fabs(c.molecules - 10000) < 10);
  /* the duplicates of all reads: 20000 reads of 10000 molecules */
  test_assert(fabs(c.molecules_seen - 8646.6) < 10);
  /* no duplicates */
  complexity_estimate(1000, 1000, 24, 0.0, &c);
  test_assert(c.saturated_unknown);
  /* 10% errors: p = 0.9^24 = 0.08, no estimate */
  complexity_estimate(20000, 19000, 24, 0.10, &c);
  test_assert(!c.estimated && c.high_error_rate);
  test_assert(fabs(c.error_free - pow(0.9, 24)) < 1e-9);
  /* ... but too few reads */
  complexity_estimate(50, 50, 24, 0.10, &c);
  test_assert(!c.estimated && !c.high_error_rate);
  /* 3%: estimated */
  test_assert(pow(0.97, 24) >= COMPLEXITY_MIN_ERROR_FREE);
}
