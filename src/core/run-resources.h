/* (c) 2026 Giorgio Gonnella */

/*
   The resources a run uses, for the dashboard and the report: elapsed
   time, CPU time, memory, and how many threads are busy.

   The CPU time is that of the whole process (all its threads, user and
   system). The load is the CPU time divided by the elapsed time, i.e.
   the number of threads busy on average: between two measures while the
   run is in progress (the current load), from the start at its end. The
   memory is the resident memory of the process now (Linux: from
   /proc/self/statm; 0 where it cannot be read) and its peak: the greater
   of the peak the system reports (getrusage, which it may update with a
   delay) and the memory measured now or at an earlier measure, so that
   the peak is never below the memory now.
*/

#ifndef RUN_RESOURCES_H
#define RUN_RESOURCES_H

#include <stdbool.h>
#include <stddef.h>

/* a measure of the resources of the run */
typedef struct {
  size_t threads;        /* the threads of the trimming (-j) */
  double elapsed;        /* seconds since the start of the run */
  double cpu;            /* CPU seconds used since the start */
  double load;           /* threads busy on average, see above */
  size_t memory;         /* resident memory now (bytes; 0: unknown) */
  size_t memory_peak;    /* peak resident memory (bytes) */
} run_resources_t;

/* the clock of a run: its start, the previous measure, and the greatest
   memory measured */
typedef struct {
  double start, last_elapsed, last_cpu;
  size_t memory_peak;
} run_clock_t;

/* Start the clock of a run (now) */
void run_clock_start(run_clock_t *clock);

/* Measure the resources of the run of <threads> threads; the load is
   that since the previous measure (since the start for the first one,
   or if <whole_run>) */
void run_resources_measure(run_clock_t *clock, size_t threads,
    bool whole_run, run_resources_t *out);

/* "1:02:03" (hours:minutes:seconds) or "2:03" (minutes:seconds) */
typedef struct { char s[32]; } run_duration_t;
run_duration_t run_duration_text(double seconds);

/* A memory size in bytes, formatted as e.g. "850 MB" or "1.24 GB". */
typedef struct { char s[24]; } run_memory_text_t;
run_memory_text_t run_memory_text(size_t bytes);

/* Unit tests */
void run_resources_test(void);

#endif
