/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include "html-page.h"
#include "html-report-markup.h"
#include "html-report-run.h"
#include "html-template.h"
#include "numutils.h"
#include "platform-check.h"
#include "segments.h"

/* the rows of the table of the settings (html-report-*.html), embedded by
   the Makefile */
static const char *const run_tmpl =
#include "html-report-run-html.inc"
  ;
static const char *const options_tmpl =
#include "html-report-options-html.inc"
  ;

/* the errors in a primer: the option, or decided automatically */
static void write_max_errors(html_writer_t *w, const report_t *r) {
  const report_setting_t *s = r->max_errors;
  double value = s != NULL ? s->value : r->options->max_errors;
  markup_open_row(w);
  markup_td(w, "errors in a primer");
  markup_open(w, "td", "seq");
  if (s != NULL && s->pending)
    html_text(w, "automatic, at the first input file");
  else {
    if (value < 1)
      html_textf(w, "%g%% of its length", 100 * value);
    else
      html_textf(w, "%g", value);
    if (s != NULL && s->automatic)
      markup_textf(w, "span", "note", " (automatic: error rate %s in the "
                   "first %s reads of %s)",
                   fraction_percent_text(s->measured).s,
                   html_num(s->n_reads).s, s->file);
  }
  html_close(w);
  html_close(w);
}

/* the primer scheme, if chosen automatically (-s auto) */
static void write_scheme_choice(html_writer_t *w, const report_t *r) {
  const scheme_choice_t *c = r->scheme_choice;
  if (c == NULL)
    return;
  markup_open_row(w);
  markup_td(w, "scheme choice");
  markup_open(w, "td", "seq");
  if (c->pending)
    html_text(w, "automatic, at the first input file");
  else {
    html_textf(w, "automatic: %s, in the first %s reads of %s",
               scheme_choice_outcome_text(c->outcome),
               html_num(c->n_reads).s, c->file);
    for (size_t e = 0; e < 2; e++) {
      char text[256];
      scheme_choice_end_text(c, e, text, sizeof(text));
      html_elem(w, "br");
      markup_textf(w, "span", "note", "%s", text);
    }
  }
  html_close(w);
  html_close(w);
}

/* the check of the platform against the reads (platform-check.h) */
static void write_platform_check(html_writer_t *w, const report_t *r) {
  const platform_check_t *c = r->platform_check;
  if (c == NULL)
    return;
  markup_open_row(w);
  markup_td(w, "platform check");
  markup_open(w, "td", "seq");
  if (c->pending)
    html_text(w, c->automatic ? "automatic, at the first input file" :
              "at the first input file");
  else {
    char evidence[512];
    platform_check_evidence(c, evidence, sizeof(evidence));
    if (c->reads == PLATFORM_READS_FEW)
      html_textf(w, "not done: %s", evidence);
    else if (c->mismatch)
      html_textf(w, "the reads may be of the other platform: %s",
                 evidence);
    else if (c->automatic)
      html_textf(w, "platform read from the reads: %s", evidence);
    else
      html_textf(w, "the reads agree: %s", evidence);
    if (!c->automatic && !c->stop)
      markup_textf(w, "span", "note", " (not stopping on a contradiction: "
                   "--no-platform-check)");
  }
  html_close(w);
  html_close(w);
}

/* where the reads look like the other platform than the one used, a
   warning above the settings */
static void write_platform_warning(html_writer_t *w, const report_t *r) {
  const platform_check_t *c = r->platform_check;
  if (c == NULL || c->pending || !c->mismatch)
    return;
  char evidence[512];
  platform_check_evidence(c, evidence, sizeof(evidence));
  const char *used = c->illumina ? "Illumina" : "nanopore";
  const char *other = c->illumina ? "nanopore" : "Illumina";
  markup_open(w, "p", "platformwarn");
  html_tag_text(w, "b", "Platform");
  html_textf(w, " \xc2\xb7 these reads may be %s reads, but the %s "
             "settings were used: %s. If they are %s reads, run again "
             "with -p %s.", other, used, evidence, other,
             c->illumina ? "nanopore" : "illumina");
  html_close(w);
}

/* all the settings, as a table */
static void write_settings(html_writer_t *w, const report_t *r) {
  const options_t *o = r->options;
  char n_refs[32], window[32], min_length[32], probe[32], probe_errors[32],
       terminal_window[32], threads[32], segment[32];
  if (o->report_segment != 0)
    segments_name(o->report_segment, virus_set_main(o->virus), segment,
                  sizeof(segment));
  snprintf(n_refs, sizeof(n_refs), "%zu", r->references == NULL ? 0 :
           r->references->n_records);
  snprintf(window, sizeof(window), "%ld", o->window);
  snprintf(min_length, sizeof(min_length), "%zu", o->min_length);
  snprintf(probe, sizeof(probe), "%zu", o->probe_len);
  snprintf(probe_errors, sizeof(probe_errors), "%g%% of its length",
           100 * o->probe_error_rate);
  snprintf(terminal_window, sizeof(terminal_window), "%ld",
           o->terminal_window);
  snprintf(threads, sizeof(threads), "%zu", o->threads);
  html_open(w, "table");
  html_template(w, run_tmpl, "version", r->version,
                "command", r->command_line, "scheme", r->scheme->name,
                "policy", o->policy == TRIM_PRIMER ?
                  "primer (every primer-derived base removed)" :
                  "terminus (the conserved termini kept)",
                "references", r->references != NULL ? o->references : NULL,
                "n_references", n_refs,
                "no_references", r->references == NULL ? "" : NULL);
  write_scheme_choice(w, r);
  write_max_errors(w, r);
  write_platform_check(w, r);
  html_template(w, options_tmpl, "platform",
                o->illumina ? "illumina" : "nanopore", "window", window,
                "min_length", min_length,
                "probe", r->references != NULL ? probe : NULL,
                "probe_errors", probe_errors,
                "terminal_window", terminal_window,
                "search", o->exhaustive_search ? ", exhaustive search" : "",
                "segment", o->segment != 0 ? segment : NULL,
                "report_segment", o->report_segment != 0 ? segment : NULL,
                "threads", threads);
  html_close(w);
}

/* the main settings, on one line */
static void write_line(html_writer_t *w, const report_t *r) {
  markup_open(w, "p", "runline");
  html_tag_text(w, "b", r->options->illumina ? "Illumina" : "Nanopore");
  html_textf(w, " \xc2\xb7 primers %s \xc2\xb7 %s", r->scheme->name,
             r->options->policy == TRIM_PRIMER ? "primers removed whole" :
             "conserved termini kept");
  if (r->references != NULL)
    html_textf(w, " \xc2\xb7 %zu references", r->references->n_records);
  else
    html_text(w, " \xc2\xb7 no references");
  if (r->options->report_segment != 0) {
    char segment[32];
    segments_name(r->options->report_segment,
                  virus_set_main(r->options->virus), segment,
                  sizeof(segment));
    html_textf(w, " \xc2\xb7 %s only", segment);
  }
  html_textf(w, " \xc2\xb7 flucut %s", r->version);
  html_close(w);
  /* the sequencing run, if the reads of all the samples are of one */
  const read_metadata_t *m = counter_metadata(r->total);
  const char *cell = read_metadata_text(m, READ_METADATA_FLOW_CELL);
  const char *model = read_metadata_text(m, READ_METADATA_MODEL);
  const char *name = read_metadata_text(m, READ_METADATA_EXPERIMENT);
  bool one_cell = cell != NULL && !m->several[READ_METADATA_FLOW_CELL];
  bool one_model = model != NULL && !m->several[READ_METADATA_MODEL];
  bool one_name = name != NULL && !m->several[READ_METADATA_EXPERIMENT];
  if (one_cell || one_model || one_name) {
    markup_open(w, "p", "runline");
    if (one_name)
      html_textf(w, "%s", name);
    if (one_cell)
      html_textf(w, "%sflow cell %s", one_name ? " \xc2\xb7 " : "", cell);
    if (one_model)
      html_textf(w, "%s%s", one_name || one_cell ? " \xc2\xb7 " : "",
                 model);
    html_close(w);
  }
}

void run_write(html_writer_t *w, const report_t *r) {
  assert(r != NULL);
  if (r->running) {
    markup_open(w, "p", "running");
    html_tag_text(w, "b", "Run in progress");
    html_textf(w, " \xc2\xb7 updated after each %s%s",
               r->options->live ? "input file" : "sample",
               r->refresh > 0 ? " \xc2\xb7 this page reloads itself" : "");
    html_close(w);
  } else if (r->dashboard) {
    markup_open(w, "p", "finished");
    html_tag_text(w, "b", "Run finished");
    html_text(w, " \xc2\xb7 all samples on one page, for printing: ");
    html_open(w, "a");
    html_attr(w, "href", "%s", r->report_file);
    html_text(w, r->report_file);
    html_close(w);
    html_close(w);
  }
  write_line(w, r);
  write_platform_warning(w, r);
  markup_details_start(w, "Settings");
  write_settings(w, r);
  markup_details_end(w);
}
