/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "alloc.h"
#include "progress.h"
#include "strutils.h"
#include "test.h"

/* width of the bar, in characters */
#define PROGRESS_BAR_WIDTH 20

struct progress_t {
  FILE *file;
  char *description;
  size_t total, count;
  double start, last_refresh;
  size_t last_len;     /* length of the last line, for overwriting it */
  bool finished;
};

/* seconds from an arbitrary start (the wall clock, with a resolution of a
   second, if the monotonic clock is not available) */
static double now(void) {
  struct timespec ts;
  if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0)
    return (double)time(NULL);
  return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

/* format a time interval as [h:]mm:ss (truncated), as tqdm; not
   format_duration() of strutils.h, whose format is [D d ]HH:MM:SS
   (rounded), e.g. 00:01:05 instead of 01:05 */
static void format_time(char *buf, size_t size, double seconds) {
  unsigned long s = seconds < 0 ? 0 : (unsigned long)seconds;
  if (s >= 3600)
    snprintf(buf, size, "%lu:%02lu:%02lu", s / 3600, (s / 60) % 60, s % 60);
  else
    snprintf(buf, size, "%02lu:%02lu", s / 60, s % 60);
}

/* blocks for the partial character of the bar (1/8 to 7/8) */
static const char *partial_blocks[] = {"", "▏", "▎", "▍",
  "▌", "▋", "▊", "▉"};
#define FULL_BLOCK "█"

static void refresh(progress_t *self) {
  char line[512], elapsed[32], remaining[32];
  double t = now() - self->start;
  double rate = t > 0 ? (double)self->count / t : 0;
  format_time(elapsed, sizeof(elapsed), t);
  size_t len;
  if (self->total == 0)
    len = (size_t)snprintf(line, sizeof(line), "%s: %zuit [%s, %.2fit/s]",
                           self->description, self->count, elapsed, rate);
  else {
    double fraction = (double)self->count / (double)self->total;
    if (fraction > 1)
      fraction = 1;
    if (rate > 0)
      format_time(remaining, sizeof(remaining),
                  (double)(self->total - (self->count < self->total ?
                           self->count : self->total)) / rate);
    else
      strcpy(remaining, "?");
    char bar[PROGRESS_BAR_WIDTH * 3 + 1] = "";
    size_t eighths = (size_t)(fraction * PROGRESS_BAR_WIDTH * 8);
    size_t full = eighths / 8, chars = 0;
    for (size_t i = 0; i < full; i++, chars++)
      strcat(bar, FULL_BLOCK);
    if (full < PROGRESS_BAR_WIDTH) {
      strcat(bar, partial_blocks[eighths % 8]);
      if (eighths % 8 > 0)
        chars++;
    }
    for (; chars < PROGRESS_BAR_WIDTH; chars++)
      strcat(bar, " ");
    len = (size_t)snprintf(line, sizeof(line),
                           "%s: %3.0f%%|%s| %zu/%zu [%s<%s, %.2fit/s]",
                           self->description, fraction * 100, bar,
                           self->count, self->total, elapsed, remaining,
                           rate);
  }
  if (len >= sizeof(line))
    len = sizeof(line) - 1;
  fprintf(self->file, "\r%s", line);
  /* clear the rest of a previous longer line */
  for (size_t i = len; i < self->last_len; i++)
    fputc(' ', self->file);
  self->last_len = len;
  fflush(self->file);
  self->last_refresh = now();
}

progress_t* progress_new(FILE *file, const char *description, size_t total) {
  assert(file != NULL);
  assert(description != NULL);
  progress_t *self;
  XMALLOC(self, sizeof(progress_t));
  self->file = file;
  self->description = xstrdup(description);
  self->total = total;
  self->count = 0;
  self->start = now();
  self->last_len = 0;
  self->finished = false;
  refresh(self);
  return self;
}

void progress_update(progress_t *self, size_t n) {
  assert(self != NULL);
  self->count += n;
  if (now() - self->last_refresh >= PROGRESS_MIN_INTERVAL)
    refresh(self);
}

void progress_finish(progress_t *self) {
  assert(self != NULL);
  if (self->finished)
    return;
  refresh(self);
  fputc('\n', self->file);
  fflush(self->file);
  self->finished = true;
}

size_t progress_count(const progress_t *self) {
  assert(self != NULL);
  return self->count;
}

void progress_delete(progress_t *self) {
  if (self == NULL)
    return;
  free(self->description);
  free(self);
}

void progress_test(void) {
  char *buf;
  size_t size;
  FILE *f = open_memstream(&buf, &size);
  progress_t *p = progress_new(f, "Counting", 0);
  progress_update(p, 5);
  progress_update(p, 7);
  test_eq_size(progress_count(p), 12);
  progress_finish(p);
  progress_delete(p);
  fclose(f);
  /* the final line */
  char *last = strrchr(buf, '\r');
  test_assert(last != NULL && strncmp(last, "\rCounting: 12it [00:00, ", 24)
              == 0);
  test_assert(buf[size - 1] == '\n');
  free(buf);
  f = open_memstream(&buf, &size);
  p = progress_new(f, "Bar", 8);
  progress_update(p, 4);
  progress_finish(p);
  progress_delete(p);
  fclose(f);
  last = strrchr(buf, '\r');
  /* half of the bar (10 full blocks of 3 bytes, 10 spaces) */
  test_assert(strncmp(last, "\rBar:  50%|", 11) == 0);
  test_assert(strstr(last, "| 4/8 [00:00<") != NULL);
  free(buf);
}
