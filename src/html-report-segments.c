/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include "html-page.h"
#include "html-report-markup.h"
#include "html-report-plots.h"
#include "html-report-sample.h"
#include "html-report-segments.h"
#include "numutils.h"
#include "segments.h"

/* the main internal priming site of a segment is marked on its coverage
   from this share of the depth from internally primed reads */
#define SITE_MARK_SHARE 0.05

/* a number of the strip of a segment, with its name */
static void number(html_writer_t *w, const char *name, const char *value) {
  html_open(w, "div");
  html_tag_text(w, "dt", name);
  html_tag_text(w, "dd", value);
  html_close(w);
}

/* the numbers of the reference measured as <q> */
static void write_numbers(html_writer_t *w, const qc_segment_t *q) {
  char depth[40], full[40], name[40];
  snprintf(depth, sizeof(depth), "%s\xc3\x97",
           html_num(q->median_depth).s);
  snprintf(full, sizeof(full), "%s\xc3\x97",
           html_num(q->median_full_depth).s);
  snprintf(name, sizeof(name), "at \xe2\x89\xa5%d\xc3\x97", QC_MIN_DEPTH);
  markup_open(w, "dl", "kv");
  number(w, "median depth", depth);
  number(w, "of full templates", full);
  number(w, name, fraction_percent_text(q->well_covered).s);
  number(w, "covered", fraction_percent_text(q->covered).s);
  markup_textf(w, "div", "kh", "share of the depth");
  number(w, "internally primed", fraction_percent_text(q->internal).s);
  number(w, "internal deletion", fraction_percent_text(q->deletion).s);
  number(w, "minus strand", fraction_percent_text(q->minus).s);
  html_close(w);
}

/* the name of a strip, linked to the part <sample> of the report (the
   index of a sample; -1: no link) */
static void write_name(html_writer_t *w, const char *name, long sample) {
  html_open(w, "b");
  if (sample >= 0) {
    html_open(w, "a");
    html_attr(w, "href", "#s%ld", sample);
    html_text(w, name);
    html_close(w);
  } else
    html_text(w, name);
  html_close(w);
}

/* the state of the segment measured as <q>, as a badge */
static void segment_badge(html_writer_t *w, const qc_segment_t *q) {
  markup_badge(w, q->status == QC_SEGMENT_COMPLETE ? MARKUP_GOOD :
               q->status == QC_SEGMENT_PARTIAL ? MARKUP_WARNING :
               MARKUP_CRITICAL, "%s",
               q->status == QC_SEGMENT_COMPLETE ? "complete" :
               q->status == QC_SEGMENT_PARTIAL ? "partial" : "low depth");
}

/* the name of a strip with the state of the segment (<badge>) */
static void write_title(html_writer_t *w, const qc_segment_t *q,
    const char *name, long sample, bool badge) {
  markup_open(w, "div", "sname");
  write_name(w, name, sample);
  if (badge)
    segment_badge(w, q);
  html_close(w);
}

/* the strip of the reference measured as <q>, named <name> (a gene, the
   reference itself, or a sample: see write_name()); <badge>: with the
   state of the segment; <subtype>: the subtype of the sample, shown
   under the name (NULL: not shown). In the dashboard the strip is
   folded: only its name and its main numbers are shown until it is
   opened; in the report (which is printed) it is always open */
static void write_strip(html_writer_t *w, const report_t *r,
    const counter_t *c, const qc_segment_t *q, const char *name,
    long sample, bool badge, const char *subtype) {
  const fasta_record_t *ref = r->references->records + q->ref;
  bool fold = r->dashboard;
  if (fold) {
    html_open(w, "details");
    html_attr(w, "class", "segd");
    html_open(w, "summary");
    write_title(w, q, name, sample, badge);
    if (subtype != NULL && subtype[0] != '\0')
      markup_textf(w, "span", "sbrief", "%s", subtype);
    markup_textf(w, "span", "sbrief", "median depth %s\xc3\x97 (full "
                 "templates %s\xc3\x97) \xc2\xb7 %s at \xe2\x89\xa5%d\xc3\x97 "
                 "\xc2\xb7 internally primed %s", html_num(q->median_depth).s,
                 html_num(q->median_full_depth).s,
                 fraction_percent_text(q->well_covered).s, QC_MIN_DEPTH,
                 fraction_percent_text(q->internal).s);
    html_close(w);
  }
  markup_open(w, "div", "seg");
  markup_open(w, "div", "sinfo");
  if (!fold)
    write_title(w, q, name, sample, badge);
  if (!fold && subtype != NULL && subtype[0] != '\0')
    markup_textf(w, "div", "ssub", "%s", subtype);
  markup_open(w, "div", "note");
  if (badge)
    html_textf(w, "%s \xc2\xb7 ", ref->name);
  if (q->subtype[0] != '\0')
    html_textf(w, "%s \xc2\xb7 ", q->subtype);
  html_textf(w, "%s b", html_num(q->ref_len).s);
  html_close(w);
  write_numbers(w, q);
  html_close(w);
  markup_open(w, "div", "scov");
  plots_coverage(w, c, q->ref, q->ref_len, QC_MIN_DEPTH,
                 q->internal >= SITE_MARK_SHARE);
  html_close(w);
  markup_open(w, "div", "slen");
  if (q->complete > 0) {
    plots_lengths(w, counter_lengths(c), q->ref, q->ref_len);
    markup_note(w, "%s complete amplicons, %s short",
                html_num(q->complete).s,
                percent_text(q->short_complete, q->complete).s);
  } else
    markup_note(w, "No complete amplicons.");
  html_close(w);
  html_close(w);
  if (fold)
    html_close(w);
}

/* the strip of a segment without located reads, or (<waiting>) of a
   sample not processed yet */
static void write_missing(html_writer_t *w, const char *name, long sample,
    bool waiting) {
  markup_open(w, "div", "seg missing");
  markup_open(w, "div", "sinfo");
  markup_open(w, "div", "sname");
  write_name(w, name, sample);
  if (waiting)
    markup_badge(w, MARKUP_NEUTRAL, "waiting");
  else
    markup_badge(w, MARKUP_CRITICAL, "not found");
  html_close(w);
  html_close(w);
  html_close(w);
}

/* the titles and the legends of the plots, above the strips */
static void write_head(html_writer_t *w) {
  markup_open(w, "div", "seghead");
  html_tag_text(w, "span", "");
  html_open(w, "div");
  html_tag_text(w, "b", "Coverage by read origin");
  plots_coverage_legend(w, QC_MIN_DEPTH);
  html_close(w);
  html_open(w, "div");
  html_tag_text(w, "b", "Complete amplicons by length (bases)");
  markup_open(w, "div", "lg");
  markup_swatch(w, "c-full", "full length");
  markup_swatch(w, "c-short", "short");
  html_close(w);
  html_close(w);
  html_close(w);
}

void segments_write(html_writer_t *w, const report_t *r, const counter_t *c,
    const qc_summary_t *qc) {
  assert(r->references != NULL);
  int shown = r->options->report_segment;
  markup_open(w, "section", "card");
  html_tag_text(w, "h3", "Segments");
  write_head(w);
  for (int s = 1; s <= SEGMENTS_N; s++) {
    if ((shown != 0 && shown != s) || !qc_expects(qc, s))
      continue;
    if (qc->segments[s].found)
      write_strip(w, r, c, qc->segments + s, segments_gene(s), -1, true,
                  NULL);
    else
      write_missing(w, segments_gene(s), -1, false);
  }
  /* the references whose description names no segment */
  for (size_t g = 0; g < r->references->n_records && shown == 0; g++) {
    const fasta_record_t *ref = r->references->records + g;
    if (segments_of_description(ref->description) != 0 ||
        !counter_has_coverage(c, g))
      continue;
    qc_segment_t q;
    qc_reference(r->references, c, g, &q);
    write_strip(w, r, c, &q, ref->name, -1, false, NULL);
  }
  /* (the last group of the lengths: the amplicons not located) */
  lengths_summary_t not_located;
  lengths_summarize(counter_lengths(c), r->references->n_records, 0,
                    &not_located);
  if (not_located.total > 0)
    markup_note(w, "%s complete amplicons were not located in the "
                "references.", html_num(not_located.total).s);
  html_close(w);
}

/* the segment whose cleavage site is read from the reads */
#define CLEAVAGE_SEGMENT 4

/* the cell of the HA cleavage site: the residues the reads carry and
   how far they can be trusted (cleavage-site.h) */
static void cleavage_cell(html_writer_t *w, const report_t *r,
    const counter_t *c) {
  cleavage_result_t site;
  report_cleavage(r, c, &site);
  html_open(w, "td");
  if (!site.found)
    html_text(w, "\xe2\x80\x93");
  else {
    html_open(w, "span");
    html_attr(w, "class", "motif-cell");
    html_textf(w, "%s", site.aa);
    html_close(w);
    markup_badge(w,
        site.verdict == CLEAVAGE_FIRM ? MARKUP_GOOD :
        site.verdict == CLEAVAGE_TOO_FEW_READS ? MARKUP_CRITICAL :
        MARKUP_WARNING, "%s", cleavage_verdict_name(site.verdict));
  }
  html_close(w);
}

/* a row of the table of one segment in all samples: the sample <i> */
static void across_row(html_writer_t *w, const report_t *r, size_t i,
    int segment) {
  const report_section_t *g = r->sections + i;
  markup_open_row(w);
  html_open(w, "td");
  html_open(w, "a");
  html_attr(w, "href", "#s%zu", i);
  html_text(w, g->name);
  html_close(w);
  html_close(w);
  qc_summary_t qc;
  if (g->counts != NULL)
    qc_summary(r->references, g->counts, r->options->segment, &qc);
  if (g->counts == NULL || !qc.segments[segment].found) {
    html_open(w, "td");
    html_attr(w, "class", "wait");
    html_attr(w, "colspan", segment == CLEAVAGE_SEGMENT ? "11" : "10");
    html_text(w, g->counts == NULL ? "waiting" : "not found");
    html_close(w);
    html_close(w);
    return;
  }
  const qc_segment_t *q = qc.segments + segment;
  html_open(w, "td");
  segment_badge(w, q);
  html_close(w);
  markup_td(w, "%s", sample_subtype(&qc).s);
  markup_td(w, "%s\xc3\x97", html_num(q->median_depth).s);
  markup_td(w, "%s\xc3\x97", html_num(q->median_full_depth).s);
  markup_td(w, "%s", fraction_percent_text(q->well_covered).s);
  markup_td(w, "%s", fraction_percent_text(q->covered).s);
  markup_td(w, "%s", fraction_percent_text(q->internal).s);
  markup_td(w, "%s", fraction_percent_text(q->deletion).s);
  markup_td(w, "%s", html_num(q->complete).s);
  markup_td(w, "%s", q->complete > 0 ?
            percent_text(q->short_complete, q->complete).s : "\xe2\x80\x93");
  if (segment == CLEAVAGE_SEGMENT)
    cleavage_cell(w, r, g->counts);
  html_close(w);
}

/* the table of one segment in all samples, above their strips */
static void write_across_table(html_writer_t *w, const report_t *r,
    int segment) {
  char well[32];
  snprintf(well, sizeof(well), "at \xe2\x89\xa5%d\xc3\x97", QC_MIN_DEPTH);
  markup_open(w, "div", "scroll");
  markup_open(w, "table", "samples across");
  if (segment == CLEAVAGE_SEGMENT)
    markup_header_row(w, "sample", "state", "subtype", "median depth",
                      "of full templates", well,
                      "covered", "int. primed", "int. deletions",
                      "complete amplicons", "short", "cleavage site", NULL);
  else
    markup_header_row(w, "sample", "state", "subtype", "median depth",
                      "of full templates", well,
                      "covered", "int. primed", "int. deletions",
                      "complete amplicons", "short", NULL);
  for (size_t i = 0; i < r->n_sections; i++)
    across_row(w, r, i, segment);
  html_close(w);
  html_close(w);
  markup_note(w, "of full templates: the median depth counting only the "
              "reads whose amplicon started at a segment end. "
              "%s: share of the %s positions at a depth of %d or more; "
              "covered: at a depth of 1 or more; int. primed, int. "
              "deletions: shares of the depth from internally primed reads "
              "and from reads with an internal deletion.", well,
              segments_gene(segment), QC_MIN_DEPTH);
  if (segment == CLEAVAGE_SEGMENT)
    markup_note(w, "cleavage site: the residues before the fusion peptide "
                "as the full templates carry them, read from the reads "
                "and not from an assembly; the badge says how far they "
                "can be trusted (see the sample below).");
}

void segments_write_across(html_writer_t *w, const report_t *r,
    int segment) {
  assert(r->references != NULL);
  markup_open(w, "section", "card");
  write_across_table(w, r, segment);
  write_head(w);
  for (size_t i = 0; i < r->n_sections; i++) {
    const report_section_t *g = r->sections + i;
    if (g->counts == NULL) {
      write_missing(w, g->name, (long)i, true);
      continue;
    }
    qc_summary_t qc;
    qc_summary(r->references, g->counts, r->options->segment, &qc);
    if (qc.segments[segment].found)
      write_strip(w, r, g->counts, qc.segments + segment, g->name, (long)i,
                  true, sample_subtype(&qc).s);
    else
      write_missing(w, g->name, (long)i, false);
  }
  html_close(w);
}

void segments_write_balance(html_writer_t *w, const qc_summary_t *qc) {
  double total = 0, max = 0;
  for (int s = 1; s <= SEGMENTS_N; s++) {
    if (!qc->segments[s].found)
      continue;
    total += qc->segments[s].bases;
    if (qc->segments[s].bases > max)
      max = qc->segments[s].bases;
  }
  if (total <= 0)
    return;
  markup_open(w, "section", "card");
  html_tag_text(w, "h3", "Bases by segment");
  markup_open(w, "div", "columns");
  for (int s = 1; s <= SEGMENTS_N; s++) {
    double bases = qc->segments[s].found ? qc->segments[s].bases : 0;
    html_open(w, "div");
    markup_textf(w, "span", "cv", "%.0f%%", 100 * bases / total);
    html_open(w, "i");
    html_attr(w, "style", "height:%.0fpx", 70 * bases / max);
    html_close(w);
    html_tag_text(w, "b", segments_gene(s));
    html_close(w);
  }
  html_close(w);
  markup_note(w, "Share of the bases of the located reads.");
  html_close(w);
}

void segments_write_lengths(html_writer_t *w, const counter_t *c) {
  markup_open(w, "section", "card");
  html_tag_text(w, "h3", "Complete amplicons by length (bases)");
  plots_lengths(w, counter_lengths(c), 0, 0);
  markup_note(w, "Without references the read ends are not classified: "
              "internal priming shows here as amplicons shorter than their "
              "segment.");
  html_close(w);
}
