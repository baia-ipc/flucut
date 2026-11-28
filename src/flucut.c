/* (c) 2026 Giorgio Gonnella */

/*
   flucut: removal of the primers of the influenza A multisegment RT-PCR
   from the reads, before the assembly (see README.md).

   Each group of input files (input-groups.h: a file of reads, a
   directory of files of reads, or each barcode directory of a MinKNOW
   fastq_pass directory) gives one gzipped FASTQ file in the output
   directory, and a line of the summary table (summary.tsv); the
   reports (HTML and JSON) are written at the end, and, with several
   groups, also after each of them. The files are trimmed
   by worker threads (worker.h), whose output does not depend on their
   number. With --live, the files are processed as they are written
   (live.h).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "alloc.h"
#include "auto-settings.h"
#include "banner.h"
#include "builtin-references.h"
#include "markers.h"
#include "counter.h"
#include "fasta.h"
#include "file.h"
#include "group-state.h"
#include "input-groups.h"
#include "live.h"
#include "logger.h"
#include "options.h"
#include "platform-check.h"
#include "read-outputs.h"
#include "reference.h"
#include "results.h"
#include "run-report.h"
#include "scheme.h"
#include "scheme-choice.h"
#include "segments.h"
#include "trimmer.h"
#include "worker.h"


/*
   Whether the references look like coding sequences (CDS) without the
   untranslated regions (UTRs) of the segments: most start with a start
   codon and end with a stop codon, and none begins with the conserved
   terminus of a segment. flucut locates the read in the references by the
   sequence just after the 5' primer (the probe), which in a full segment
   is the 5' UTR; a CDS reference does not contain it, so the terminal
   reads cannot be classified (they are reported as not located), and only
   the internally primed reads, whose probe falls inside the ORF, map.
   The references must be complete segments, with the UTRs.
*/
static bool references_look_cds_only(const fasta_t *refs) {
  if (refs == NULL || refs->n_records == 0)
    return false;
  size_t cds = 0;
  for (size_t i = 0; i < refs->n_records; i++) {
    const char *s = refs->records[i].seq;
    size_t n = refs->records[i].seq_len;
    if (n < 6)
      continue;
    /* a complete segment begins with the universal terminus (AGCrAAAGCAGG);
       such a reference is never flagged */
    if (n >= 12 && (memcmp(s, "AGCAAAAGCAGG", 12) == 0 ||
                    memcmp(s, "AGCGAAAGCAGG", 12) == 0))
      return false;
    bool starts = memcmp(s, "ATG", 3) == 0;
    const char *e = s + n - 3;
    bool ends = memcmp(e, "TAA", 3) == 0 || memcmp(e, "TAG", 3) == 0 ||
                memcmp(e, "TGA", 3) == 0;
    if (starts && ends)
      cds++;
  }
  return cds * 2 >= refs->n_records;   /* at least half look like a CDS */
}


/* the outputs of the groups (see read-outputs.h); an input which would
   be overwritten by an output, or an output file of two groups (e.g.
   those by category of the group a and the output of a group a.full),
   are errors */
static read_outputs_t** group_outputs(const input_group_t *groups,
    size_t n_groups, const options_t *opts, Logger *logger) {
  read_outputs_t **outputs;
  XMALLOC(outputs, sizeof(read_outputs_t *) * (n_groups + 1));
  for (size_t i = 0; i < n_groups; i++) {
    outputs[i] = read_outputs_new(opts->outdir, groups[i].name,
                                  opts->by_category, opts->filtered);
    for (size_t j = 0; j < groups[i].n_files; j++)
      LCHECK(!read_outputs_has_path(outputs[i], groups[i].files[j]), logger,
             "%s: the input would be overwritten by the output\n",
             groups[i].files[j]);
    for (size_t j = 0; j < i; j++)
      for (size_t k = 0; k < outputs[i]->n; k++)
        LCHECK(outputs[i]->paths[k] == NULL ||
               !read_outputs_has_path(outputs[j], outputs[i]->paths[k]),
               logger, "%s: the output file of two inputs (%s and %s)\n",
               outputs[i]->paths[k], groups[j].name, groups[i].name);
  }
  return outputs;
}

static void delete_group_outputs(read_outputs_t **outputs, size_t n) {
  for (size_t i = 0; i < n; i++)
    read_outputs_delete(outputs[i]);
  free(outputs);
}

/* say which primers are searched for a sample, where the scheme holds
   the primers of both types and its references are of one of them */
static void log_sample_virus(virus_set_t types, const char *name,
    Logger *logger) {
  if (types != VIRUS_SET_BOTH)
    lprintf(logger, LOG_INFO, "%s: influenza %s (from its references): "
            "its primers alone are searched\n", name,
            virus_set_name(types));
}

/* log the scheme (verbose) */
static void log_scheme(const scheme_t *scheme, Logger *logger) {
  lprintf(logger, LOG_INFO, "Scheme %s: left", scheme->name);
  for (size_t i = 0; i < scheme->n_left; i++)
    lnprintf(logger, LOG_INFO, " %s (%s)", scheme->left[i].name,
             scheme->left[i].seq);
  lnprintf(logger, LOG_INFO, "; right");
  for (size_t i = 0; i < scheme->n_right; i++)
    lnprintf(logger, LOG_INFO, " %s (%zu bases)", scheme->right[i].name,
             scheme->right[i].seq_len);
  lnprintf(logger, LOG_INFO, "\n");
}

/* log the input groups (verbose) */
static void log_groups(const input_group_t *groups, size_t n_groups,
    Logger *logger) {
  lprintf(logger, LOG_INFO, "Input groups: %zu\n", n_groups);
  for (size_t i = 0; i < n_groups; i++)
    lprintf(logger, LOG_INFO, "  %s: %zu file%s\n", groups[i].name,
            groups[i].n_files, groups[i].n_files == 1 ? "" : "s");
}

/* trim the files of the group <g> (the <i>-th of <n_groups>) to <out>,
   adding the counts to <counts> */
static void trim_group(const options_t *opts, const trimmer_config_t *config,
    worker_t *workers, const input_group_t *g, size_t i, size_t n_groups,
    read_outputs_t *out, counter_t *counts, Logger *logger) {
  /* (one line per group; the files of a group of several: verbose) */
  char of[64] = "";
  if (n_groups > 1)
    snprintf(of, sizeof(of), " (%zu of %zu)", i + 1, n_groups);
  if (g->n_files == 1)
    lprintf(logger, LOG_SUCCESS, "Trimming %s%s...\n", g->files[0], of);
  else
    lprintf(logger, LOG_SUCCESS, "Trimming %s, %zu files%s...\n", g->name,
            g->n_files, of);
  read_outputs_open(out, false, logger);
  /* each barcode is a different sample: classify it only against the
     reference of each segment closest to it (chosen from its first
     file), so the other references of a segment, e.g. of other
     subtypes, do not match it spuriously */
  if (g->n_files > 0) {
    /* the type of the sample, from the references its reads were
       located in: its reads are trimmed with the primers of that type
       alone */
    virus_set_t types;
    subtype_reads_t subtype_reads;
    bool *allowed = auto_settings_choose_references(config, g->files[0],
        &types, &subtype_reads, logger);
    counter_set_subtype_reads(counts, &subtype_reads);
    log_sample_virus(types, g->name, logger);
    for (size_t t = 0; t < opts->threads; t++) {
      worker_restrict_references(workers + t, allowed);
      worker_restrict_virus(workers + t, types);
    }
    free(allowed);
  }
  for (size_t j = 0; j < g->n_files; j++) {
    if (g->n_files > 1)
      lprintf(logger, LOG_PHASE, "Trimming %s...\n", g->files[j]);
    worker_process_file(workers, opts->threads, g->files[j], out, counts,
                        logger);
  }
  read_outputs_close(out, logger);
}

/* with --resume: whether the group <g> was completely trimmed by an
   earlier run in the same output directory (its saved list of files is
   the list of its files, and its outputs exist); if so, its saved counts
   are added to <counts> */
static bool group_already_trimmed(const options_t *opts,
    const input_group_t *g, const read_outputs_t *out, counter_t *counts) {
  size_t n;
  char **files = group_state_load_files(opts->outdir, g->name, &n);
  bool same = files != NULL && n == g->n_files;
  for (size_t j = 0; same && j < n; j++)
    same = strcmp(files[j], g->files[j]) == 0;
  group_state_free_files(files, n);
  for (size_t k = 0; same && k < out->n; k++)
    same = out->paths[k] == NULL || access(out->paths[k], F_OK) == 0;
  return same && group_state_load_counts(opts->outdir, g->name, counts);
}

/* trim the files of the input groups (not the live mode); with several
   groups, the tables and reports are also written before the first group
   and after each one (the dashboard, see run-report.h), so
   that a long run can be followed */
static void trim_groups(options_t *opts, trimmer_config_t *config,
    platform_check_t *platform, scheme_choice_t *choice,
    auto_settings_t *autos, scheme_t *scheme, const fasta_t *refs,
    run_report_t *report, Logger *logger) {
  size_t n_groups;
  input_group_t *groups = input_groups_new(opts->inputs, opts->n_inputs,
                                           &n_groups, false, true, logger);
  read_outputs_t **outputs = group_outputs(groups, n_groups, opts, logger);
  log_groups(groups, n_groups, logger);
  /* the platform (whose window the scheme choice uses), -s auto, then
     -e auto (which measures the primers of the scheme): from the first
     file */
  for (size_t i = 0; i < n_groups && platform->pending; i++)
    if (groups[i].n_files > 0)
      platform_check_decide(platform, opts, config, autos,
                            groups[i].files[0], logger);
  for (size_t i = 0; i < n_groups && choice->pending; i++)
    if (groups[i].n_files > 0)
      scheme_choice_decide(choice, scheme, opts->virus, opts->window,
                           groups[i].files[0], logger);
  for (size_t i = 0; i < n_groups && auto_settings_pending(autos); i++)
    if (groups[i].n_files > 0)
      auto_settings_decide(autos, config, groups[i].files[0], logger);
  worker_t *workers;
  XMALLOC(workers, sizeof(worker_t) * opts->threads);
  for (size_t t = 0; t < opts->threads; t++)
    worker_init(workers + t, config, opts->compress_level);
  /* the counts of each group (NULL: not processed yet), and of all reads */
  counter_t *total = counter_new(scheme, refs);
  results_group_t *results;
  XMALLOC(results, sizeof(results_group_t) * (n_groups + 1));
  for (size_t i = 0; i < n_groups; i++)
    results[i] = (results_group_t){groups[i].name, outputs[i]->label, NULL,
                                   groups[i].n_files};
  for (size_t i = 0; i < n_groups; i++) {
    if (n_groups > 1)
      results_write(opts->outdir, results, n_groups, total, refs, report,
                    true, logger);
    if (n_groups > 1 && i == 0)
      run_report_announce_dashboard(report, logger);
    counter_t *counts = counter_new(scheme, refs);
    if (opts->resume &&
        group_already_trimmed(opts, groups + i, outputs[i], counts)) {
      lprintf(logger, LOG_SUCCESS, "%s: trimmed by an earlier run, "
              "skipped (--resume)\n", groups[i].name);
    } else {
      /* (the counts of a group whose loading failed half-way are
         discarded; its state is removed, so that it is only saved again
         once the group is complete) */
      counter_delete(counts);
      counts = counter_new(scheme, refs);
      group_state_remove(opts->outdir, groups[i].name, logger);
      trim_group(opts, config, workers, groups + i, i, n_groups, outputs[i],
                 counts, logger);
      group_state_save(opts->outdir, groups[i].name, counts,
                       groups[i].files, groups[i].n_files, logger);
    }
    counter_merge(total, counts);
    results[i].counts = counts;
  }
  results_write(opts->outdir, results, n_groups, total, refs, report,
                false, logger);
  /* (the states of the groups are kept: --resume on a complete run only
     writes the reports again) */
  for (size_t i = 0; i < n_groups; i++)
    counter_delete((counter_t *)results[i].counts);
  free(results);
  counter_delete(total);
  for (size_t t = 0; t < opts->threads; t++)
    worker_free(workers + t);
  free(workers);
  delete_group_outputs(outputs, n_groups);
  input_groups_delete(groups, n_groups);
}

/* warn if the classification (--segment) or the report
   (--report-segment) are limited to a segment which no reference is
   recognised as (see segments.h) */
static void check_segment(const options_t *opts, const fasta_t *refs,
    Logger *logger) {
  if (opts->report_segment == 0 || refs == NULL)
    return;
  for (size_t i = 0; i < refs->n_records; i++)
    if (report_shows_reference(refs, opts->report_segment, i))
      return;
  lprintf(logger, LOG_WARN, "%s: no reference is recognised as segment %d "
          "(%s) from its description; %s\n",
          opts->segment != 0 ? "--segment" : "--report-segment",
          opts->report_segment, segments_gene(opts->report_segment),
          opts->segment != 0 ? "no read end will be located" :
          "the per-reference parts of the report will be empty");
}

int main(int argc, char *argv[]) {
  options_t opts;
  Logger *logger;
  options_parse(&opts, argc, argv, &logger);
  if (opts.show_banner)
    banner_print(stderr);
  /* -s auto: the default of --virus until the scheme is chosen from the
     first reads (scheme-choice.h) */
  scheme_choice_t choice;
  scheme_choice_init(&choice, &opts);
  scheme_t *scheme = scheme_new_builtin(scheme_choice_initial(&opts),
                                        logger);
  if (choice.automatic)
    lprintf(logger, LOG_INFO, "Scheme: chosen from the first reads (-s "
            "auto)\n");
  else
    log_scheme(scheme, logger);
  /* the references: -r off none; -r <file> the user's; the default the
     built-in set of the types expected (--virus), downloaded on first
     use (builtin-references.h) */
  const char *ref_path = opts.references_off ? NULL :
      opts.references != NULL ? opts.references :
      builtin_references_path(opts.virus, logger);
  fasta_t *refs = ref_path == NULL ? NULL : fasta_read(ref_path, logger);
  /* the reports show the actual path used (the built-in cache, or the
     user's file) */
  opts.references = ref_path;
  if (refs == NULL)
    opts.policy = TRIM_PRIMER;   /* no references: every primer removed */
  reference_index_t *index = NULL;
  if (refs != NULL) {
    lprintf(logger, LOG_INFO, "%zu reference sequences read from %s\n",
            refs->n_records, ref_path);
    if (references_look_cds_only(refs))
      lprintf(logger, LOG_WARN,
          "the references look like coding sequences without the UTRs "
          "(they start at a start codon and end at a stop codon); the 5' "
          "origin of the reads cannot be classified (they will be reported "
          "as not located) and only internal priming inside the ORF is "
          "seen. Use complete segment references, with the UTRs "
          "(scripts/fetch_references.py).\n");
    index = reference_index_new(refs);
  }
  check_segment(&opts, refs, logger);
  /* the cleavage site of the haemagglutinin is read from the reads
     (cleavage-site.h); it needs an HA reference which contains the
     fusion peptide */
  cleavage_anchors_t *cleavage = cleavage_anchors_new(refs, index, logger);
  if (cleavage != NULL)
    lprintf(logger, LOG_INFO,
            "the cleavage site of the HA will be read from the reads "
            "(%zu HA references with the fusion peptide)\n",
            cleavage_anchors_n(cleavage));
  /* the markers are read from the reads beside the conserved peptides
     which anchor them (markers.h); a reference gives a marker only
     where its anchor is in it exactly once */
  marker_anchors_t *markers = marker_anchors_new(refs, index,
                                                 opts.probe_error_rate,
                                                 logger);
  if (markers != NULL)
    lprintf(logger, LOG_INFO,
            "the markers will be read from the reads (%zu anchors in the "
            "references)\n", marker_anchors_n(markers));
  /* the platform: checked against the first reads, or read from them
     (-p auto); until then, the settings of the platform given, or of
     nanopore (platform-check.h) */
  platform_check_t platform;
  platform_check_init(&platform, &opts);
  auto_settings_t autos;
  auto_settings_init(&autos, &opts);
  trimmer_config_t config = {.scheme = scheme,
      .max_errors = autos.max_errors.value, .policy = opts.policy,
      .window = opts.window, .references = refs, .index = index,
      .exhaustive = opts.exhaustive_search,
      .seeds_only = opts.illumina && !opts.exhaustive_search,
      .segment = opts.segment,
      .probe_len = opts.probe_len,
      .probe_error_rate = opts.probe_error_rate,
      .terminal_window = opts.terminal_window,
      .adapters = opts.adapters, .min_length = opts.min_length,
      .illumina = opts.illumina, .cleavage = cleavage,
      .markers = markers,
      .tags = opts.tags, .chimeras = opts.chimeras,
      .deletions = opts.deletions, .filtered = opts.filtered,
      .by_category = opts.by_category};
  run_report_t *report = run_report_new(&opts, VERSION, scheme, refs, argc,
                                        argv, logger);
  run_report_set_max_errors(report, &autos.max_errors);
  run_report_set_platform_check(report, &platform);
  if (choice.automatic)
    run_report_set_scheme_choice(report, &choice);
  if (opts.live)
    live_run(&opts, &config, &platform, &choice, &autos, scheme, refs,
             report, logger);
  else
    trim_groups(&opts, &config, &platform, &choice, &autos, scheme, refs,
                report, logger);
  run_report_delete(report);
  cleavage_anchors_delete(cleavage);
  marker_anchors_delete(markers);
  auto_settings_free(&autos);
  platform_check_free(&platform);
  scheme_choice_free(&choice);
  reference_index_delete(index);
  fasta_delete(refs);
  scheme_delete(scheme);
  delete_logger(logger);
  return EXIT_SUCCESS;
}
