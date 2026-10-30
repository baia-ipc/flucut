/* (c) 2026 Giorgio Gonnella */

/*
   Content of the reports of a run (HTML report, JSON summary), as in
   Lumipore: the counts of all reads and of each group of input files,
   the information needed to describe them, and the values computed from
   them for both reports (the internal priming sites, the short complete
   amplicons).
*/

#ifndef REPORT_H
#define REPORT_H

#include <stdbool.h>
#include <stddef.h>
#include "counter.h"
#include "fasta.h"
#include "lengths.h"
#include "options.h"
#include "scheme.h"

/* suffix of the temporary file of a report, renamed when complete */
#define REPORT_TMP_SFX ".tmp"

/* the counts of a group of input files */
typedef struct {
  const char *name;
  const counter_t *counts;  /* NULL: not processed (yet) */
  size_t n_files;
} report_section_t;

/* a setting which can be decided automatically (e.g. the errors in a
   primer, -e auto): the value; if decided automatically (and not
   <pending>), the measurement (e.g. the error rate) in the first
   <n_reads> reads of <file> */
typedef struct {
  bool automatic, pending;
  double value, measured;
  size_t n_reads;
  const char *file;
} report_setting_t;

typedef struct {
  const char *version;
  const char *command_line;
  const options_t *options;
  const scheme_t *scheme;
  const fasta_t *references;  /* NULL: none */
  const report_section_t *sections;
  size_t n_sections;
  const counter_t *total;
  /* the errors in a primer, if decided automatically (NULL: the value of
     the options) */
  const report_setting_t *max_errors;
  /* a run in progress: reload the page every <refresh> seconds (0: no);
     <dashboard>: a dashboard (one screen), instead of a report (not
     written yet: the report is written) */
  long refresh;
  bool running;
  bool dashboard;
} report_t;

/* the name of an end (COUNTER_LEFT, COUNTER_RIGHT) and that of the
   origin of a read end (end_origin_t) */
const char* report_end_name(size_t end);
const char* report_origin_name(end_origin_t origin);

/* the read ends at positions at most this far apart are one internal
   priming site (the probes of the reads of a site map a few bases
   apart, because of the errors) */
#define REPORT_SITE_GAP 5

/* an internal priming site: the read ends of an end (COUNTER_LEFT or
   COUNTER_RIGHT) classified as internal, mapped to the positions <from>
   to <to> (1-based, see counter_sites()) of the reference <ref>, most of
   them to <pos> */
typedef struct {
  size_t ref, end;
  size_t pos, from, to;
  size_t reads;
} report_site_t;

/*
   The internal priming sites of the counts <c> (for the references
   <refs>, NULL: none), by decreasing number of reads (then by reference,
   end and position); their number is stored in <n>; the array is
   allocated (to be freed by the caller; NULL if there are none)
*/
report_site_t* report_sites(const fasta_t *refs, const counter_t *c,
    size_t *n);

/* the size of the label of a reference (see segments_label()) */
#define REPORT_LABEL_SIZE 256

/* the complete amplicons of a group of the lengths (see
   counter_lengths()): with references, of the reference <group> (or not
   located, if <group> is the number of references), with the amplicons
   shorter than LENGTHS_SHORT_FRACTION of the reference (none for the
   amplicons not located) */
typedef struct {
  const char *name;        /* the reference, "not_located" or "all" */
  char label[REPORT_LABEL_SIZE]; /* the name, with the segment of the
                                    reference if known (segments.h) */
  size_t ref_len;          /* 0: no reference */
  lengths_summary_t summary;
} report_lengths_t;

/* The complete amplicons of the group <group> of the lengths of <c> */
void report_lengths(const fasta_t *refs, const counter_t *c, size_t group,
    report_lengths_t *out);

/* The complete amplicons shorter than LENGTHS_SHORT_FRACTION of their
   reference, in all groups of <c> (0 without references) */
size_t report_short_complete(const fasta_t *refs, const counter_t *c);

/* Unit tests */
void report_test(void);

#endif
