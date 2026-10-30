/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdlib.h>
#include "alloc.h"
#include "file.h"
#include "html-report.h"
#include "json-report.h"
#include "run-report.h"
#include "sam-header.h"

struct run_report_t {
  const options_t *opts;
  const char *version;
  const scheme_t *scheme;
  const fasta_t *references;
  const report_setting_t *max_errors;
  char *command_line, *html_path, *json_path;
};

run_report_t* run_report_new(const options_t *opts, const char *version,
    const scheme_t *scheme, const fasta_t *references, int argc,
    char *argv[], Logger *logger) {
  assert(opts != NULL);
  assert(scheme != NULL);
  run_report_t *self;
  XCALLOC(self, 1, sizeof(run_report_t));
  self->opts = opts;
  self->version = version;
  self->scheme = scheme;
  self->references = references;
  /* (as in the @PG line of a BAM file) */
  self->command_line = sam_header_command_line(argc, argv);
  xmakedirs(opts->outdir, logger);
  self->html_path = path_join(opts->outdir, RUN_REPORT_HTML_FILE);
  self->json_path = path_join(opts->outdir, RUN_REPORT_JSON_FILE);
  return self;
}

void run_report_set_max_errors(run_report_t *self,
    const report_setting_t *setting) {
  assert(self != NULL);
  self->max_errors = setting;
}

void run_report_write(const run_report_t *self,
    const report_section_t *sections, size_t n, const counter_t *total,
    bool running, Logger *logger) {
  assert(self != NULL);
  report_t r = {
    .version = self->version, .command_line = self->command_line,
    .options = self->opts, .scheme = self->scheme,
    .references = self->references, .sections = sections,
    .n_sections = n, .total = total, .max_errors = self->max_errors,
    .refresh = running ? RUN_REPORT_REFRESH : 0, .running = running,
    .dashboard = false};
  html_report_write(self->html_path, &r, logger);
  json_report_write(self->json_path, &r, logger);
}

void run_report_delete(run_report_t *self) {
  if (self == NULL)
    return;
  free(self->command_line);
  free(self->html_path);
  free(self->json_path);
  free(self);
}
