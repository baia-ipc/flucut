/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <errno.h>
#include <glob.h>
#include <signal.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include "khash.h"
#include "alloc.h"
#include "counter.h"
#include "file.h"
#include "input-groups.h"
#include "live.h"
#include "reads-reader.h"
#include "report.h"
#include "results.h"
#include "strutils.h"
#include "worker.h"

/* the longest single sleep while waiting for the next scan (seconds),
   so that a signal (live_stop) is noticed quickly */
#define LIVE_SLEEP_STEP 0.1

KHASH_MAP_INIT_STR(strsize, long)

/* set by the signal handler (SIGINT, SIGTERM); file-level, as a signal
   handler can only communicate through such a variable (the only state
   of the module) */
static volatile sig_atomic_t live_stop = 0;

static void live_signal(int sig) {
  (void)sig;
  live_stop = 1;
}

/* state of a group in the live mode */
typedef struct {
  char *name;
  char *output;            /* the path of its output */
  char **done;             /* the files processed, in processing order */
  size_t n_done, done_cap;
  khash_t(strsize) *done_set;   /* file sizes (-1: unknown/warned) */
  counter_t *counts;
} live_group_t;

static bool is_done(const live_group_t *st, const char *path) {
  return kh_get(strsize, st->done_set, path) != kh_end(st->done_set);
}

static long path_size(const char *path) {
  struct stat st;
  return stat(path, &st) == 0 ? (long)st.st_size : -1L;
}

/* warn (once) if a processed file was modified after its processing */
static void check_size(live_group_t *st, const char *path, Logger *logger) {
  khiter_t k = kh_get(strsize, st->done_set, path);
  assert(k != kh_end(st->done_set));
  long size = kh_value(st->done_set, k);
  if (size >= 0 && path_size(path) != size) {
    lprintf(logger, LOG_WARN, "%s was modified after being processed: "
            "the new reads are not processed (increase --min-age)\n",
            path);
    kh_value(st->done_set, k) = -1;
  }
}

static void add_done(live_group_t *st, const char *path, long size) {
  if (st->n_done == st->done_cap) {
    st->done_cap = st->done_cap == 0 ? 64 : st->done_cap * 2;
    XREALLOC(st->done, sizeof(char *) * st->done_cap);
  }
  char *copy = xstrdup(path);
  st->done[st->n_done++] = copy;
  int ret;
  khiter_t k = kh_put(strsize, st->done_set, copy, &ret);
  if (ret == -1)
    fatal_printf("Failed allocating the processed files\n");
  kh_value(st->done_set, k) = size;
}

/* the path of a state file of a group (to be freed by the caller) */
static char* state_path(const char *outdir, const char *group,
    const char *sfx) {
  char *name = strcat_alloc(group, sfx);
  char *path = path_join(outdir, name);
  free(name);
  return path;
}

/* write a state file of a group, the counts or the list of processed
   files: to a temporary file, then renamed, so that an interruption
   while writing it does not leave it incomplete */
static void write_state_file(const options_t *opts, const live_group_t *st,
    bool counts, Logger *logger) {
  char *path = state_path(opts->outdir, st->name,
                          counts ? LIVE_COUNTS_SFX : LIVE_FILES_SFX);
  FILE *f = xfopen(path, REPORT_TMP_SFX, "w", "state file", logger);
  if (counts)
    counter_save(st->counts, f);
  else
    for (size_t i = 0; i < st->n_done; i++)
      fprintf(f, "%s\n", st->done[i]);
  xfclose(&f, logger);
  xfrename(path, REPORT_TMP_SFX, path, NULL, logger);
  free(path);
}

/*
   save the counts and the list of processed files of a group (for
   --resume); the list is replaced last, so that a file is only listed as
   processed once its counts are saved: after an interruption between the
   two, the file is processed again on resuming (as after an interruption
   while processing it; then its counts are added again, as its reads are
   appended again to the output), instead of being skipped with its
   counts missing
*/
static void save_state(const options_t *opts, const live_group_t *st,
    Logger *logger) {
  write_state_file(opts, st, true, logger);
  write_state_file(opts, st, false, logger);
}

/* with --resume: read the list of processed files and the counts of a
   group, if available */
static void load_state(const options_t *opts, live_group_t *st,
    Logger *logger) {
  char *path = state_path(opts->outdir, st->name, LIVE_FILES_SFX);
  FILE *f = fopen(path, "r");
  free(path);
  if (f == NULL)
    return;
  char *line = NULL;
  size_t cap = 0;
  ssize_t len;
  while ((len = getline(&line, &cap, f)) > 0) {
    if (line[len - 1] == '\n')
      line[--len] = '\0';
    if (len > 0 && !is_done(st, line))
      add_done(st, line, -1);
  }
  free(line);
  fclose(f);
  path = state_path(opts->outdir, st->name, LIVE_COUNTS_SFX);
  f = fopen(path, "r");
  free(path);
  if (f == NULL || !counter_load(st->counts, f))
    lprintf(logger, LOG_WARN, "%s: the counts of the files processed "
            "before resuming could not be read\n", st->name);
  if (f != NULL)
    fclose(f);
  lprintf(logger, LOG_INFO, "%s: %zu file%s processed before resuming\n",
          st->name, st->n_done, st->n_done == 1 ? "" : "s");
}

/* remove the file <path>, if any */
static void remove_file(const char *path, Logger *logger) {
  LCHECK(unlink(path) == 0 || errno == ENOENT, logger,
         "Deleting file '%s' failed, %s\n", path, strerror(errno));
}

/* remove the output of a group and its state files (when it is processed
   from the beginning, as the output is then appended to) */
static void clear_group(const options_t *opts, const live_group_t *st,
    Logger *logger) {
  remove_file(st->output, logger);
  const char *sfx[] = {LIVE_FILES_SFX, LIVE_COUNTS_SFX};
  for (size_t i = 0; i < 2; i++) {
    char *path = state_path(opts->outdir, st->name, sfx[i]);
    remove_file(path, logger);
    free(path);
  }
}

/* the live mode: the arguments of live_run() and its state */
typedef struct {
  const options_t *opts;
  trimmer_config_t *config;
  auto_settings_t *autos;
  const scheme_t *scheme;
  const fasta_t *refs;
  run_report_t *report;
  Logger *logger;
  worker_t *workers;       /* NULL: not yet created (see get_workers()) */
  live_group_t *groups;    /* the groups found until now */
  size_t n_groups;
  counter_t *file_counts;  /* the counts of the file being processed */
  bool ended;              /* the end of the run was found */
  long min_age;            /* minimum age of the files of the scan */
  size_t n_incomplete;     /* files incomplete at the end of the run */
  time_t last_new;         /* when the last file was processed */
} live_t;

/* the state of the group <name>, created (and, with --resume, loaded) if
   new */
static live_group_t* get_group(live_t *w, const char *name) {
  const options_t *opts = w->opts;
  for (size_t i = 0; i < w->n_groups; i++)
    if (strcmp(w->groups[i].name, name) == 0)
      return w->groups + i;
  XREALLOC(w->groups, sizeof(live_group_t) * (w->n_groups + 1));
  live_group_t *st = w->groups + w->n_groups++;
  st->name = xstrdup(name);
  st->output = input_group_output_path(opts->outdir, name);
  st->done = NULL;
  st->n_done = 0;
  st->done_cap = 0;
  st->done_set = kh_init(strsize);
  st->counts = counter_new(w->scheme, w->refs);
  if (opts->resume)
    load_state(opts, st, w->logger);
  if (st->n_done == 0)
    clear_group(opts, st, w->logger);
  return st;
}

/*
   a file is ready to be processed if it has not been modified for the
   minimum age, and its content is complete (for gzipped files: the gzip
   stream is not truncated; FASTQ: the last line is complete; BAM: the
   file ends with the BGZF end-of-file marker)
*/
static bool file_ready(const char *path, long min_age) {
  struct stat st;
  /* (an empty file is being created, e.g. by a copy) */
  if (stat(path, &st) != 0 || st.st_size == 0 ||
      difftime(time(NULL), st.st_mtime) < min_age)
    return false;
  return reads_reader_file_complete(path);
}

/* whether a file of the end of the run (LIVE_END_FILES) is in <dir> */
static bool end_file_in(const char *dir) {
  char *pattern = strcat_alloc(dir, "/" LIVE_END_FILES);
  glob_t g;
  bool found = glob(pattern, 0, NULL, &g) == 0 && g.gl_pathc > 0;
  globfree(&g);
  free(pattern);
  return found;
}

/* whether the run has ended: a file LIVE_END_FILES in an input directory
   or in its parent (e.g. the run directory of fastq_pass) */
static bool run_ended(const options_t *opts) {
  for (size_t i = 0; i < opts->n_inputs; i++) {
    char *parent = strcat_alloc(opts->inputs[i], "/..");
    bool found = end_file_in(opts->inputs[i]) || end_file_in(parent);
    free(parent);
    if (found)
      return true;
  }
  return false;
}

/* write the tables and reports of the groups; <running>: before the end
   of the mode */
static void write_results(const live_t *w, bool running) {
  results_group_t *results;
  XMALLOC(results, sizeof(results_group_t) * (w->n_groups + 1));
  counter_t *total = counter_new(w->scheme, w->refs);
  for (size_t i = 0; i < w->n_groups; i++) {
    const live_group_t *st = w->groups + i;
    results[i] = (results_group_t){st->name, st->output, st->counts,
                                   st->n_done};
    counter_merge(total, st->counts);
  }
  results_write(w->opts->outdir, results, w->n_groups, total, w->refs,
                w->report, running, w->logger);
  counter_delete(total);
  free(results);
}

/* stop at SIGINT and SIGTERM (see live_stop) */
static void set_signals(Logger *logger) {
  struct sigaction sa;
  memset(&sa, 0, sizeof(sa));
  sa.sa_handler = live_signal;
  /* SA_RESETHAND: the handler is reset to the default action after the
     first signal, so that a second Ctrl-C kills the process, e.g. if
     the processing of a large file takes too long; the first one only
     sets live_stop (the outputs and reports are then completed) */
  /* (the cast: SA_RESETHAND is unsigned, e.g. 0x80000000) */
  sa.sa_flags = (int)SA_RESETHAND;
  LCHECK(sigaction(SIGINT, &sa, NULL) == 0 &&
         sigaction(SIGTERM, &sa, NULL) == 0, logger, "Setting the handler "
         "of SIGINT and SIGTERM failed, %s\n", strerror(errno));
}

/* -e auto: decided from the first file ready (also if already processed,
   so that a resumed run uses the same) */
static void decide_auto(live_t *w, const char *path) {
  if (!auto_settings_pending(w->autos) || !file_ready(path, w->min_age))
    return;
  auto_settings_decide(w->autos, w->config, path, w->logger);
}

/* the workers, created at the first file processed (after -e auto) */
static worker_t* get_workers(live_t *w) {
  if (w->workers == NULL) {
    XMALLOC(w->workers, sizeof(worker_t) * w->opts->threads);
    for (size_t t = 0; t < w->opts->threads; t++)
      worker_init(w->workers + t, w->config);
  }
  return w->workers;
}

/* process a new file <path> of size <size> of a group, appending to its
   output, and save its state, the tables and the reports */
static void process_file(live_t *w, live_group_t *st, const char *path,
    long size) {
  const options_t *opts = w->opts;
  counter_clear(w->file_counts);
  FILE *out = xfopen(st->output, NULL, st->n_done > 0 ? "ab" : "wb",
                     "trimmed reads", w->logger);
  worker_process_file(get_workers(w), opts->threads, path, out, st->output,
                      w->file_counts, w->logger);
  xfclose(&out, w->logger);
  counter_merge(st->counts, w->file_counts);
  add_done(st, path, size);
  save_state(opts, st, w->logger);
  const char *base = strrchr(path, '/');
  char *name = strcat_alloc(st->name, ": ");
  char *label = strcat_alloc(name, base == NULL ? path : base + 1);
  results_log_counts(label, counter_reads(w->file_counts), w->refs != NULL,
                     w->logger);
  free(label);
  free(name);
  write_results(w, true);
  w->last_new = time(NULL);
}

/* the files of a group found in a scan: the new ones which are ready are
   processed */
static void scan_group(live_t *w, const input_group_t *group) {
  live_group_t *st = get_group(w, group->name);
  for (size_t i = 0; i < group->n_files && !live_stop; i++) {
    const char *path = group->files[i];
    decide_auto(w, path);
    if (is_done(st, path)) {
      check_size(st, path, w->logger);
      continue;
    }
    /* the size is checked to be the same before and after file_ready():
       otherwise data were appended while checking that the content is
       complete; the size is that of the content processed (see
       check_size()) */
    long size = path_size(path);
    if (!file_ready(path, w->min_age) || path_size(path) != size) {
      w->n_incomplete += w->ended;
      continue;
    }
    process_file(w, st, path, size);
  }
}

/* scan the input directories; false if the mode must end */
static bool scan(live_t *w) {
  const options_t *opts = w->opts;
  /* the end of the run is checked before the scan, so that the files
     written before it are processed; then they are complete: no minimum
     age */
  w->ended = run_ended(opts);
  w->min_age = w->ended ? 0 : opts->live_min_age;
  size_t n_groups;
  input_group_t *groups = input_groups_new(opts->inputs, opts->n_inputs,
      &n_groups, true, false, w->logger);
  for (size_t g = 0; g < n_groups && !live_stop; g++)
    scan_group(w, groups + g);
  input_groups_delete(groups, n_groups);
  return !w->ended && !live_stop && (opts->live_timeout == 0 ||
      difftime(time(NULL), w->last_new) < opts->live_timeout);
}

/* log why the mode ended */
static void log_end(const live_t *w) {
  if (live_stop)
    lprintf(w->logger, LOG_SUCCESS, "Interrupted: stopped watching\n");
  else if (w->ended) {
    lprintf(w->logger, LOG_SUCCESS, "End of the run (" LIVE_END_FILES "): "
            "stopped watching\n");
    if (w->n_incomplete > 0)
      lprintf(w->logger, LOG_WARN, "%zu file%s incomplete at the end of the "
              "run, not processed\n", w->n_incomplete,
              w->n_incomplete > 1 ? "s were" : " was");
  } else
    lprintf(w->logger, LOG_SUCCESS, "No new files for %ld s: stopped "
            "watching\n", w->opts->live_timeout);
}

static void delete_groups(live_group_t *groups, size_t n_groups) {
  for (size_t i = 0; i < n_groups; i++) {
    for (size_t j = 0; j < groups[i].n_done; j++)
      free(groups[i].done[j]);
    free(groups[i].done);
    kh_destroy(strsize, groups[i].done_set);
    counter_delete(groups[i].counts);
    free(groups[i].output);
    free(groups[i].name);
  }
  free(groups);
}

/* sleep <seconds>, in steps of at most LIVE_SLEEP_STEP, until a signal
   sets live_stop */
static void live_sleep(double seconds) {
  while (seconds > 0 && !live_stop) {
    double step = seconds < LIVE_SLEEP_STEP ? seconds : LIVE_SLEEP_STEP;
    struct timespec ts = {.tv_sec = (time_t)step,
                          .tv_nsec = (long)((step - (time_t)step) * 1e9)};
    nanosleep(&ts, NULL);
    seconds -= step;
  }
}

void live_run(const options_t *opts, trimmer_config_t *config,
    auto_settings_t *autos, const scheme_t *scheme,
    const fasta_t *references, run_report_t *report, Logger *logger) {
  assert(opts != NULL);
  assert(config != NULL);
  assert(autos != NULL);
  set_signals(logger);
  xmakedirs(opts->outdir, logger);
  if (opts->live_timeout > 0)
    lprintf(logger, LOG_PHASE, "Watching the input directories (stop: end "
            "of the run (" LIVE_END_FILES "), %ld s without new files, "
            "or Ctrl-C)\n", opts->live_timeout);
  else
    lprintf(logger, LOG_PHASE, "Watching the input directories (stop: end "
            "of the run (" LIVE_END_FILES "), or Ctrl-C)\n");
  live_t w = {.opts = opts, .config = config, .autos = autos,
      .scheme = scheme, .refs = references, .report = report,
      .logger = logger, .workers = NULL, .groups = NULL, .n_groups = 0,
      .file_counts = counter_new(scheme, references), .ended = false,
      .min_age = opts->live_min_age, .n_incomplete = 0,
      .last_new = time(NULL)};
  /* the reports of a run in progress, before the first file */
  write_results(&w, true);
  double interval = opts->live_timeout > 0 &&
                    opts->live_timeout < opts->live_scan_interval ?
                    (double)opts->live_timeout : opts->live_scan_interval;
  while (!live_stop && scan(&w))
    live_sleep(interval);
  log_end(&w);
  counter_delete(w.file_counts);
  write_results(&w, false);
  if (w.workers != NULL) {
    for (size_t t = 0; t < opts->threads; t++)
      worker_free(w.workers + t);
    free(w.workers);
  }
  delete_groups(w.groups, w.n_groups);
}
