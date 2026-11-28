/* (c) 2026 Giorgio Gonnella */

/*
   The state of the analyzer of the reads (amplicon.h), shared by the
   modules which carry out the analysis, each one a step of it:

   - amplicon.c: the analyzer, and the analysis of a read in both
     orientations, of which the better one is kept;
   - amplicon-ends.c: the search of the primers at the ends of the read,
     and of the primer matches inside it;
   - amplicon-cut.c: the classification of the ends against the
     references, and the part of the read to keep;
   - amplicon-adapters.c: the adapters at the ends without a primer, and
     the second analysis of the part between them;
   - amplicon-locate.c: the placement of a read of which no end was
     located, the ends located from inside the read, and the junction
     of an internal deletion.

   Only these modules include this header; the rest of the program uses
   amplicon.h.
*/

#ifndef AMPLICON_ANALYZER_H
#define AMPLICON_ANALYZER_H

#include <stdbool.h>
#include <stddef.h>
#include "adapters.h"
#include "amplicon.h"
#include "internal-matches.h"
#include "primer-search.h"
#include "reference.h"
#include "scheme.h"
#include "virus.h"

struct amplicon_analyzer_t {
  const scheme_t *scheme;
  /* the errors allowed in a primer (see amplicon_analyzer_new()) */
  double max_errors;
  trim_policy_t policy;
  /* the primer of an end must lie within this many bases of that end of
     the read (see amplicon.h) */
  long window;
  /* the length of the longest primer of the scheme */
  size_t max_primer_len;
  /* the errors allowed in each primer of each end */
  long *left_max_errors, *right_max_errors;
  /* the costs of the missing outer bases of each primer of each end,
     max_primer_len + 1 for each primer (see
     amplicon_missing_costs()) */
  long *left_skip, *left_edge, *right_skip, *right_edge;
  primer_search_t *search;
  /* the search of the primers inside the read */
  internal_matches_t *internal;
  /* the types of influenza whose primers are searched (see
     amplicon_analyzer_restrict_virus) */
  virus_set_t types;
  /* the reverse complement of the read */
  char *rc;
  size_t rc_cap;
  /* the references, if any (see amplicon_analyzer_use_references()) */
  reference_set_t *refs;
  /* the segment of each reference (0: its description names none),
     read once from the descriptions: a segment-specific primer keeps
     its terminus only at an end located in a segment it binds (see
     primers.h) */
  int *ref_segment;
  long terminal_window;
  /* the probe of the right end (reference_set_probe_len() + 1 bytes) */
  char *probe;
  /* the adapters, if any (see amplicon_analyzer_use_adapters()) */
  adapters_t *adapters;
  /* the template of a read with an internal deletion, in the
     orientation of its reference (see amplicon_find_junction()) */
  char *template;
  size_t template_cap;
  /* the bases between an end and the probe which located it from inside
     (see amplicon-locate.c) */
  char *skipped;
  size_t skipped_cap;
  amplicon_t result;
};

/* The analysis of the read <seq> (length <len>), without the adapters,
   into <out>; <start_window>, <end_window>: the windows at the start and
   at the end of the read (see amplicon.c) */
void amplicon_analyze(amplicon_analyzer_t *self, const char *seq,
    size_t len, long start_window, long end_window, amplicon_t *out);

#endif
