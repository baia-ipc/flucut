/* (c) 2026 Giorgio Gonnella */

#include <stdbool.h>
#include <stdio.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>
#include "run-resources.h"
#include "test.h"

/* seconds of the monotonic clock */
static double now(void) {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (double)t.tv_sec + (double)t.tv_nsec / 1e9;
}

static double seconds(struct timeval t) {
  return (double)t.tv_sec + (double)t.tv_usec / 1e6;
}

/* the resident memory of the process now (bytes; 0: unknown) */
static size_t resident_memory(void) {
  FILE *f = fopen("/proc/self/statm", "r");
  if (f == NULL)
    return 0;
  unsigned long size = 0, resident = 0;
  bool ok = fscanf(f, "%lu %lu", &size, &resident) == 2;
  fclose(f);
  long page = sysconf(_SC_PAGESIZE);
  return ok && page > 0 ? (size_t)resident * (size_t)page : 0;
}

void run_clock_start(run_clock_t *clock) {
  clock->start = now();
  clock->last_elapsed = 0;
  clock->last_cpu = 0;
}

void run_resources_measure(run_clock_t *clock, size_t threads,
    bool whole_run, run_resources_t *out) {
  struct rusage u;
  getrusage(RUSAGE_SELF, &u);
  out->threads = threads;
  out->elapsed = now() - clock->start;
  out->cpu = seconds(u.ru_utime) + seconds(u.ru_stime);
  double from_elapsed = whole_run ? 0 : clock->last_elapsed,
         from_cpu = whole_run ? 0 : clock->last_cpu;
  double wall = out->elapsed - from_elapsed;
  out->load = wall > 0 ? (out->cpu - from_cpu) / wall : 0;
  /* (ru_maxrss: kilobytes on Linux, bytes on macOS) */
#ifdef __APPLE__
  out->memory_peak = (size_t)u.ru_maxrss;
#else
  out->memory_peak = (size_t)u.ru_maxrss * 1024;
#endif
  out->memory = resident_memory();
  clock->last_elapsed = out->elapsed;
  clock->last_cpu = out->cpu;
}

run_duration_t run_duration_text(double seconds) {
  run_duration_t text;
  long s = seconds > 0 ? (long)(seconds + 0.5) : 0;
  if (s >= 3600)
    snprintf(text.s, sizeof(text.s), "%ld:%02ld:%02ld", s / 3600,
             s / 60 % 60, s % 60);
  else
    snprintf(text.s, sizeof(text.s), "%ld:%02ld", s / 60, s % 60);
  return text;
}

void run_resources_test(void) {
  test_eq_str(run_duration_text(0).s, "0:00");
  test_eq_str(run_duration_text(59.4).s, "0:59");
  test_eq_str(run_duration_text(125).s, "2:05");
  test_eq_str(run_duration_text(3723).s, "1:02:03");
  test_eq_str(run_duration_text(-3).s, "0:00");
  /* a measure: the elapsed and CPU times grow, the memory is known */
  run_clock_t clock;
  run_clock_start(&clock);
  volatile double x = 0;
  for (long i = 0; i < 2000000; i++)
    x += (double)i;
  run_resources_t r;
  run_resources_measure(&clock, 4, false, &r);
  test_eq_size(r.threads, (size_t)4);
  test_assert(r.elapsed > 0);
  test_assert(r.cpu > 0);
  test_assert(r.load >= 0);
  test_assert(r.memory_peak > 0);
  test_assert(clock.last_elapsed == r.elapsed);
}
