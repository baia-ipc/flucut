/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "html-page.h"
#include "html-report.h"
#include "html-report-reads.h"
#include "html-report-sample.h"
#include "html-report-segments.h"
#include "html-report-subtypes.h"
#include "numutils.h"
#include "qc-flags.h"
#include "segments.h"

/* the internal priming sites shown before the reader asks for more */
#define SAMPLE_TOP_SITES 5
/* ... those with at least this number of read ends */
#define SAMPLE_MIN_SITE_READS 10

void sample_badge(html_writer_t *w, const report_t *r,
    const qc_summary_t *qc, bool brief) {
  /* (neutral: the segments cannot be judged, or not yet) */
  bool growing = r->running && r->options->live;
  markup_state_t state =
    qc->status == QC_SAMPLE_COMPLETE ? MARKUP_GOOD :
    growing || qc->status == QC_SAMPLE_FEW_READS ||
    qc->status == QC_SAMPLE_NOT_LOCATED ? MARKUP_NEUTRAL :
    qc->status == QC_SAMPLE_PARTIAL ? MARKUP_WARNING : MARKUP_CRITICAL;
  markup_badge(w, state, "%zu/%zu%s", qc->n_complete, qc->n_expected,
               brief ? "" : qc->n_expected == 1 ? " segment complete" :
               " segments complete");
}

sample_subtype_t sample_subtype(const qc_summary_t *qc) {
  /* (segment 4: HA; segment 6: NA) */
  bool unsure = qc->status > QC_SAMPLE_INCOMPLETE ||
                qc->segments[4].status <= QC_SEGMENT_LOW ||
                qc->segments[6].status <= QC_SEGMENT_LOW;
  sample_subtype_t text;
  snprintf(text.s, sizeof(text.s), "%s%s", qc->subtype, unsure ? "?" : "");
  return text;
}

sample_bases_t sample_bases(size_t bases) {
  sample_bases_t text;
  double b = (double)bases;
  if (bases >= 1000000000)
    snprintf(text.s, sizeof(text.s), "%.2f Gb", b / 1e9);
  else if (bases >= 1000000)
    snprintf(text.s, sizeof(text.s), "%.1f Mb", b / 1e6);
  else if (bases >= 1000)
    snprintf(text.s, sizeof(text.s), "%.1f kb", b / 1e3);
  else
    snprintf(text.s, sizeof(text.s), "%zu b", bases);
  return text;
}

void sample_length_tile(html_writer_t *w, const counter_t *c) {
  read_stats_summary_t in;
  char sub[96];
  read_stats_summarize(counter_read_stats(c), COUNTER_STATS_INPUT, &in);
  if (in.has_quality)
    snprintf(sub, sizeof(sub), "N50 %s, median quality Q%.0f",
             html_num(in.n50).s, in.median_quality - 0.5);
  else
    snprintf(sub, sizeof(sub), "N50 %s", html_num(in.n50).s);
  markup_tile(w, "median read length", sub, "%s",
              html_num(in.median_len).s);
}

void sample_status_line(html_writer_t *w, const report_t *r,
    const counter_t *c) {
  const counter_reads_t *s = counter_reads(c);
  bool refs = r->references != NULL;
  qc_summary_t qc;
  qc_flags_t flags;
  qc_summary(r->references, c, r->options->segment, &qc);
  cleavage_result_t site;
  report_cleavage(r, c, &site);
  qc_flags(s, &qc, refs, &site, counter_subtype_reads(c), &flags);
  if (s->reads_in == 0 || (!refs && flags.n == 0))
    return;
  markup_open(w, "p", "chips");
  if (refs) {
    sample_badge(w, r, &qc, false);
  }
  for (size_t i = 0; i < flags.n; i++)
    markup_textf(w, "span", "chip", "%s", flags.flags[i].text);
  html_close(w);
}

void sample_reads_tile(html_writer_t *w, const counter_reads_t *s,
    const run_resources_t *resources) {
  markup_open(w, "div", "tile");
  markup_textf(w, "div", "v", "%s", html_num(s->reads_in).s);
  markup_textf(w, "div", "l", "reads");
  markup_textf(w, "div", "u", "%s", sample_bases(s->bases_in).s);
  size_t analyzed = counter_reads_analyzed(s);
  if (analyzed > s->reads_in)
    markup_textf(w, "div", "u", "%s more after splitting chimeric reads",
                 html_num(analyzed - s->reads_in).s);
  if (resources != NULL && resources->elapsed > 0)
    markup_textf(w, "div", "u", "%s reads/s",
                 html_num((size_t)((double)s->reads_in /
                                   resources->elapsed + 0.5)).s);
  html_close(w);
}

void sample_written_tile(html_writer_t *w, const counter_reads_t *s) {
  size_t analyzed = counter_reads_analyzed(s);
  char sub[128];
  snprintf(sub, sizeof(sub), "%s of %s%s; %s too short",
           html_num(s->reads_out).s, html_num(analyzed).s,
           s->split > 0 ? " after splitting" : "",
           html_num(analyzed > s->reads_out ? analyzed - s->reads_out :
                    0).s);
  markup_tile(w, "written", sub, "%s",
              percent_text(s->reads_out, analyzed).s);
}

/* the main numbers of the sample */
static void write_tiles(html_writer_t *w, const report_t *r,
    const counter_t *c, const qc_summary_t *qc) {
  const counter_reads_t *s = counter_reads(c);
  size_t analyzed = counter_reads_analyzed(s);
  char sub[96];
  markup_open(w, "div", "tiles");
  if (qc->subtype[0] != '\0') {
    /* (segment 4: HA; segment 6: NA) */
    const qc_segment_t *ha = qc->segments + 4, *na = qc->segments + 6;
    snprintf(sub, sizeof(sub), "HA %s, NA %s",
             ha->found ? r->references->records[ha->ref].name : "not found",
             na->found ? r->references->records[na->ref].name : "not found");
    markup_tile(w, "subtype (closest references)", sub, "%s",
                sample_subtype(qc).s);
  }
  sample_reads_tile(w, s, NULL);
  sample_length_tile(w, c);
  sample_written_tile(w, s);
  snprintf(sub, sizeof(sub), "one %s, none %s",
           percent_text(s->n_ends[1], analyzed).s,
           percent_text(s->n_ends[0], analyzed).s);
  markup_tile(w, "both primers", sub, "%s",
              percent_text(s->n_ends[2], analyzed).s);
  html_close(w);
}

/* a count with its name, in the list of the smaller counts; <total>: it
   is also given as a share of it (0: not) */
static void small_count(html_writer_t *w, const char *name, size_t n,
    size_t total) {
  html_open(w, "div");
  html_tag_text(w, "dt", name);
  html_open(w, "dd");
  html_text(w, html_num(n).s);
  if (total > 0) {
    html_text(w, " ");
    html_tag_text(w, "small", percent_text(n, total).s);
  }
  html_close(w);
  html_close(w);
}

/* the primers found at each end, by primer of the scheme, as counts */
static void primer_counts(html_writer_t *w, const report_t *r,
    const counter_t *c, size_t total) {
  for (size_t e = 0; e < COUNTER_N_ENDS; e++) {
    const primer_t *p = e == COUNTER_LEFT ? r->scheme->left :
                        r->scheme->right;
    for (size_t v = 0; v < counter_n_variants(c, e); v++) {
      char name[96];
      snprintf(name, sizeof(name), "%s primer %s", report_end_name(e),
               p[v].name);
      small_count(w, name, counter_primer(c, e, v), total);
    }
  }
}

/* how the reads divide: by primers found, by origin (with references),
   and the smaller counts */
static void write_reads(html_writer_t *w, const report_t *r,
    const counter_t *c) {
  const counter_reads_t *s = counter_reads(c);
  size_t analyzed = counter_reads_analyzed(s);
  bool refs = r->references != NULL;
  markup_open(w, "section", "card");
  html_tag_text(w, "h3", "Reads");
  const markup_part_t primers[] = {
    {"c-p2", "both", s->n_ends[2]}, {"c-p1", "one", s->n_ends[1]},
    {"c-oth", "none", s->n_ends[0]}};
  markup_bar(w, "Primers found", "reads", primers, 3);
  if (refs) {
    const markup_part_t origin[] = {
      {"c-full", "full template", s->written[READ_CATEGORY_FULL]},
      {"c-int", "internally primed", s->written[READ_CATEGORY_INTERNAL]},
      {"c-del", "internal deletion", s->written[READ_CATEGORY_DELETION]},
      {"c-unp", "not placed", s->written[READ_CATEGORY_UNPLACED]},
      {"c-oth", "other", s->written[READ_CATEGORY_OTHER]},
      {"c-nw", "not written", s->too_short}};
    markup_bar(w, "Reads by origin", "reads", origin, 6);
  } else {
    const markup_part_t written[] = {
      {"c-full", "written", s->reads_out},
      {"c-nw", "not written", s->too_short}};
    markup_bar(w, "Reads written", "reads", written, 2);
  }
  markup_open(w, "dl", "kv");
  primer_counts(w, r, c, analyzed);
  small_count(w, "minus strand", s->minus, analyzed);
  small_count(w, "split (several amplicons)", s->split, s->reads_in);
  if (r->options->cut_deletions)
    small_count(w, "cut at the junction (internal deletion)", s->cut,
                analyzed);
  small_count(w, "primer dimers", s->dimers, analyzed);
  small_count(w, "too short", s->too_short, analyzed);
  small_count(w, "primer match inside", s->with_internal, analyzed);
  if (!r->options->illumina) {
    small_count(w, "adapter at the start", s->adapter_start, analyzed);
    small_count(w, "adapter at the end", s->adapter_end, analyzed);
  }
  html_close(w);
  if (refs) {
    markup_details_start(w, "Read ends by origin");
    for (size_t e = 0; e < COUNTER_N_ENDS; e++) {
      const markup_part_t ends[] = {
        {"c-full", "terminal", s->origin[e][ORIGIN_TERMINAL]},
        {"c-int", "internal", s->origin[e][ORIGIN_INTERNAL]},
        {"c-oth", "not located", s->origin[e][ORIGIN_UNKNOWN]}};
      markup_bar(w, e == COUNTER_LEFT ? "Left ends (Uni12 primer)" :
                 "Right ends (Uni13 primer)", "ends with a primer", ends,
                 3);
    }
    markup_details_end(w);
  }
  html_close(w);
}

/* the rows of the sites <from> to <to> */
static void site_rows(html_writer_t *w, const report_t *r,
    const report_site_t *sites, size_t from, size_t to, size_t total) {
  for (size_t i = from; i < to; i++) {
    const report_site_t *x = sites + i;
    const fasta_record_t *ref = r->references->records + x->ref;
    int segment = segments_of_description(ref->description);
    markup_open_row(w);
    markup_td(w, "%s", segment != 0 ? segments_gene(segment) : ref->name);
    markup_td(w, "%s", report_end_name(x->end));
    html_open(w, "td");
    if (x->from != x->to)
      html_attr(w, "title", "%zu-%zu", x->from, x->to);
    html_textf(w, "%s", html_num(x->pos).s);
    html_close(w);
    markup_td(w, "%s", html_num(x->reads).s);
    markup_td(w, "%s", percent_text(x->reads, total).s);
    html_close(w);
  }
}

static void sites_table(html_writer_t *w, const report_t *r,
    const report_site_t *sites, size_t from, size_t to, size_t total) {
  html_open(w, "table");
  markup_header_row(w, "segment", "end", "position", "reads",
                    "% of reads", NULL);
  site_rows(w, r, sites, from, to, total);
  html_close(w);
}

/* the internal priming sites, the largest first */
static void write_sites(html_writer_t *w, const report_t *r,
    const counter_t *c) {
  markup_open(w, "section", "card");
  html_tag_text(w, "h3", "Internal priming sites");
  size_t n, total = counter_reads_analyzed(counter_reads(c));
  report_site_t *sites = report_sites(r->references, c, &n);
  n = report_select_sites(r->references, r->options->report_segment, sites,
                          n);
  if (n == 0)
    markup_note(w, "None found.");
  else {
    size_t top = 0;
    while (top < n && top < SAMPLE_TOP_SITES &&
           sites[top].reads >= SAMPLE_MIN_SITE_READS)
      top++;
    size_t shown = n < HTML_REPORT_MAX_SITES ? n : HTML_REPORT_MAX_SITES;
    if (top > 0)
      sites_table(w, r, sites, 0, top, total);
    if (shown > top) {
      char summary[64];
      snprintf(summary, sizeof(summary), "%zu smaller sites, of %zu",
               shown - top, n);
      markup_details_start(w, summary);
      sites_table(w, r, sites, top, shown, total);
      markup_details_end(w);
    }
    markup_note(w, "Position: first base after the primer.");
  }
  free(sites);
  html_close(w);
}

/*
   The cleavage site of the haemagglutinin as the reads carry it
   (cleavage-site.h): the residues, how many reads they were read from
   and how many of them agree, and the other sites found, which show a
   sample holding the reads of more than one haemagglutinin
*/
static void write_cleavage_site(html_writer_t *w, const report_t *r,
    const counter_t *c) {
  cleavage_result_t site;
  report_cleavage(r, c, &site);
  if (site.searched == 0)
    return;
  markup_open(w, "section", "card");
  html_tag_text(w, "h3", "HA cleavage site");
  if (!site.found) {
    markup_note(w, "No read reaches it (%s reads searched).",
                html_num(site.searched).s);
    html_close(w);
    return;
  }
  markup_open(w, "p", "motif");
  html_textf(w, "%s", site.aa);
  html_close(w);
  markup_badge(w,
      site.verdict == CLEAVAGE_FIRM ? MARKUP_GOOD :
      site.verdict == CLEAVAGE_TOO_FEW_READS ? MARKUP_CRITICAL :
      MARKUP_WARNING,
      "%s", cleavage_verdict_name(site.verdict));
  markup_badge(w, site.reading == CLEAVAGE_MULTIBASIC ? MARKUP_NEUTRAL :
               MARKUP_WARNING, "%s", cleavage_reading_name(site.reading));
  html_open(w, "table");
  markup_open_row(w);
  markup_td(w, "reads it was read from");
  markup_td(w, "%s of %s searched", html_num(site.reads).s,
            html_num(site.searched).s);
  html_close(w);
  markup_open_row(w);
  markup_td(w, "reads with these residues");
  markup_td(w, "%s (%s)", html_num(site.motif_reads).s,
            percent_text(site.motif_reads, site.reads).s);
  html_close(w);
  markup_open_row(w);
  markup_td(w, "sites in at least %d%% of the reads",
            (int)(CLEAVAGE_MIN_SITE_SHARE * 100.0 + 0.5));
  markup_td(w, "%s", html_num(site.supported_proteins).s);
  html_close(w);
  markup_open_row(w);
  markup_td(w, "reads with a multibasic site");
  markup_td(w, "%s", percent_text(
      (size_t)(site.multibasic_share * (double)site.reads + 0.5),
      site.reads).s);
  html_close(w);
  if (site.control_reads > 0) {
    markup_open_row(w);
    markup_td(w, "reads agreeing beside the site");
    markup_td(w, "%s of %s", percent_text(
        (size_t)(site.control_support * (double)site.control_reads + 0.5),
        site.control_reads).s, html_num(site.control_reads).s);
    html_close(w);
  }
  markup_open_row(w);
  markup_td(w, "chemistry");
  markup_td(w, "%s", cleavage_chemistry_name(site.chemistry));
  html_close(w);
  if (site.fusion_aa[0] != '\0') {
    markup_open_row(w);
    markup_td(w, "fusion peptide in the reads");
    markup_td(w, "%s (%s of %s)", site.fusion_aa,
              percent_text((size_t)(site.fusion_support *
                                    (double)site.fusion_reads + 0.5),
                           site.fusion_reads).s,
              html_num(site.fusion_reads).s);
    html_close(w);
  }
  if (site.reference_aa[0] != '\0') {
    markup_open_row(w);
    markup_td(w, "reference: site and fusion peptide");
    markup_td(w, "%s, %s", site.reference_aa,
              site.reference_fusion[0] != '\0' ? site.reference_fusion
                                               : "\xe2\x80\x93");
    html_close(w);
  }
  html_close(w);   /* table */
  if (site.n_top > 1) {
    char summary[96];
    snprintf(summary, sizeof(summary),
             "the %zu most frequent sites in these reads", site.n_top);
    markup_details_start(w, summary);
    html_open(w, "table");
    markup_header_row(w, "site", "reads", "% of reads", "reads as", NULL);
    for (size_t i = 0; i < site.n_top; i++) {
      markup_open_row(w);
      markup_td(w, "%s", site.top[i].aa);
      markup_td(w, "%s", html_num(site.top[i].reads).s);
      markup_td(w, "%s", percent_text(site.top[i].reads, site.reads).s);
      markup_td(w, "%s", cleavage_reading_name(site.top[i].reading));
      html_close(w);
    }
    html_close(w);
    markup_details_end(w);
  }
  if (site.disagrees_with_reference)
    markup_note(w, "Most of the reads do not read like the reference of "
                "this sample (%s): it holds the reads of more than one "
                "haemagglutinin, and the site above is not necessarily "
                "the one of its own virus.", site.reference_aa);
  else if (site.control_weak)
    markup_note(w, "The reads of this sample do not agree even beside the "
                "site, where nothing is difficult (%s), so this is not a "
                "site they cannot resolve but reads too noisy, or holding "
                "more than one virus, to say. Confirm it before reporting "
                "the pathotype.", percent_text(
                    (size_t)(site.control_support *
                             (double)site.control_reads + 0.5),
                    site.control_reads).s);
  else if (site.verdict == CLEAVAGE_CONFIRM)
    markup_note(w, "Too few of the reads agree for this to stand on its "
                "own; confirm it before reporting the pathotype.");
  else if (site.verdict == CLEAVAGE_TOO_FEW_READS)
    markup_note(w, "Fewer than %d reads reach the site: nothing can be "
                "said about it.", CLEAVAGE_MIN_READS);
  markup_note(w, "The %d bases before the fusion peptide, taken from "
              "each read which covers them, translated. Read from the "
              "reads, not from an assembly. The fusion peptide itself is "
              "the anchor: the reads should agree on it, and a sample "
              "whose own differs from its reference's is not the subtype "
              "of that reference.", CLEAVAGE_SITE_NT);
  html_close(w);
}

void sample_write(html_writer_t *w, const report_t *r, const counter_t *c) {
  assert(r != NULL);
  assert(c != NULL);
  const counter_reads_t *s = counter_reads(c);
  qc_summary_t qc;
  qc_summary(r->references, c, r->options->segment, &qc);
  write_tiles(w, r, c, &qc);
  if (s->reads_in == 0) {
    markup_note(w, "No reads.");
    return;
  }
  markup_open(w, "div", "cols");
  write_reads(w, r, c);
  if (r->references != NULL) {
    html_open(w, "div");
    subtypes_write(w, c);
    write_cleavage_site(w, r, c);
    write_sites(w, r, c);
    segments_write_balance(w, &qc);
    html_close(w);
  } else
    segments_write_lengths(w, c);
  html_close(w);
  if (r->references != NULL)
    segments_write(w, r, c, &qc);
  reads_write(w, c);
}
