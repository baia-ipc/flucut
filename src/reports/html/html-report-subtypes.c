/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include "html-page.h"
#include "html-report-markup.h"
#include "html-report-subtypes.h"
#include "numutils.h"
#include "subtype-reads.h"

/* the row of the HA or of the NA (<which>) */
static void write_row(html_writer_t *w, size_t which,
    const subtype_reads_evidence_t *e) {
  char letter = subtype_reads_letter(which);
  markup_open_row(w);
  markup_td(w, "%s", subtype_reads_gene(which));
  markup_td(w, "%c%d: %s of %s reads", letter, e->leading,
            html_num(e->leading_reads).s, html_num(e->reads).s);
  if (e->second == 0) {
    markup_td(w, "none");
    markup_td_none(w);
  } else {
    markup_td(w, "%c%d: %s (%s)", letter, e->second,
              html_num(e->second_reads).s,
              percent_text(e->second_reads, e->reads).s);
    markup_td(w, "%.1f", e->second_expected);
  }
  markup_td(w, "%s (%s)", html_num(e->scattered_reads).s,
            percent_text(e->scattered_reads, e->reads).s);
  html_close(w);
}

/* a segment whose reads fit no reference of the set */
static void write_distant_note(html_writer_t *w, size_t which,
    const subtype_reads_evidence_t *e) {
  markup_note(w, "No close %s reference: in %s of %s reads (%s) the "
              "probes of one read agree on no subtype. No reference of "
              "the set is close to this segment, the subtype %c%d rests "
              "on %s reads, and the read ends of the segment are "
              "located less often. A second subtype reported beside it "
              "may be that scatter and not another virus.",
              subtype_reads_gene(which),
              html_num(e->scattered_reads).s, html_num(e->reads).s,
              percent_text(e->scattered_reads, e->reads).s,
              subtype_reads_letter(which), e->leading,
              html_num(e->leading_reads).s);
}

/* what was observed above the threshold, and what can cause it */
static void write_remark(html_writer_t *w, size_t which,
    const subtype_reads_evidence_t *e) {
  markup_note(w, "Reads of a second %s subtype: %c%d, %s of %s (%s). "
              "They can come from a mixed infection, from a contamination "
              "of the sample, or from barcode crosstalk between the "
              "samples of a run: flucut cannot tell these apart.",
              subtype_reads_gene(which), subtype_reads_letter(which),
              e->second, html_num(e->second_reads).s,
              html_num(e->reads).s,
              percent_text(e->second_reads, e->reads).s);
}

void subtypes_write(html_writer_t *w, const counter_t *c) {
  assert(c != NULL);
  subtype_reads_evidence_t e[SUBTYPE_READS_N_SEGMENTS];
  bool any = false;
  for (size_t which = 0; which < SUBTYPE_READS_N_SEGMENTS; which++) {
    subtype_reads_evidence(counter_subtype_reads(c), which, e + which);
    any = any || e[which].known;
  }
  if (!any)
    return;
  markup_open(w, "section", "card");
  html_tag_text(w, "h3", "HA and NA subtypes in the reads");
  html_open(w, "table");
  markup_header_row(w, "", "subtype of most reads", "second subtype",
                    "expected by chance", "reads agreeing on no subtype",
                    NULL);
  for (size_t which = 0; which < SUBTYPE_READS_N_SEGMENTS; which++)
    if (e[which].known)
      write_row(w, which, e + which);
  html_close(w);
  for (size_t which = 0; which < SUBTYPE_READS_N_SEGMENTS; which++) {
    if (e[which].above_threshold)
      write_remark(w, which, e + which);
    if (e[which].distant_reference)
      write_distant_note(w, which, e + which);
  }
  markup_note(w, "Of the first reads of the sample, those with several "
              "probes located in the HA or NA references: a read is of a "
              "subtype when all its probes match that subtype best.");
  html_close(w);
}
