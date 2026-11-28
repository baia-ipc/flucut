/* (c) 2026 Giorgio Gonnella */

/*
   Trimming of the reads of a batch: each read is analyzed (amplicon.h),
   the part to keep is written as a FASTQ record, and the reads are
   counted. A trimmer belongs to one thread; the threads work on disjoint
   parts of a batch, each into its own output buffer and counts, which
   the caller writes and sums in the order of the input.

   The record written
   ------------------

   The read keeps its title and its orientation; the sequence and the
   qualities are cut to the kept part. Unless disabled, SAM-like tags
   (tab-separated, as in the FASTQ written by Dorado) are appended to the
   title, telling what was found. The ends are those of the amplicon:
   for a read of the minus strand, the left end of the amplicon is at
   the end of the read.

     fs:Z:+|-|.     the strand of the amplicon found (. : no primer)
     fl:Z:<name>    the primer found at the left end (Uni12 side)
     fr:Z:<name>    the primer found at the right end (Uni13 side)
     tl:i:<n>       bases removed at the start of the read
     tr:i:<n>       bases removed at the end of the read
     ni:i:<n>       primer matches inside the read (if any)
     sp:Z:<k>/<n>   the record is the part <k> of the <n> parts of a read
                    which contained several amplicons (see below)
     jf:Z:<k>/2     the record is the flank <k> of a read cut at the
                    junction of its internal deletion (see below)

   With references, for each end found:

     ol:Z:, or:Z:   origin of the left/right end: terminal, internal or
                    unknown (see amplicon.h)
     ml:Z:, mr:Z:   where its probe mapped: <reference>:<pos>:<strand>
                    (1-based position on the plus strand)

   and for each read:

     ct:Z:<name>    its category (see read-cover.h): full (a full
                    template), internal (internal priming), deletion (an
                    internal deletion) or other (no 5' primer, or not
                    located)

   Reads by category
   -----------------

   The records can be written apart by category (<by_category>, which
   needs references), one output for each: the parts of a split read
   are then written each to the output of its own category.

   Reads shorter than the minimum length after the trimming are not
   written (primer dimers, with nothing left, always).

   Reads with an internal deletion, cut at the junction
   -----------------------------------------------------

   A read with an internal deletion is two flanks of a segment joined at
   a junction (junction.h). It can be written as two records, its two
   flanks (<cut_deletions>, which needs references), so that a program
   which maps the reads to the whole segment gets two reads which each
   lie in one place, instead of one read with a gap of hundreds of
   bases. The kept part of the read is cut at the junction: the bases
   before it in the read are the flank 1, those from it on the flank 2.
   The records have the name of the read followed by _j1 and _j2, the
   tags of the read (tl and tr are the bases removed at the two ends of
   the read, ct is deletion) and the tag jf. A flank shorter than the
   minimum length is not written. A read of the category deletion whose
   junction was not located (the ends tell the deletion, the bases
   between them do not fit the two flanks) is written whole.

   Reads with several amplicons
   ----------------------------

   A read which contains several amplicons is split at their junctions
   (chimera.h), unless disabled: each part is analyzed, trimmed, counted
   and written as a read of its own. The records of the parts have the
   name of the read followed by _<k> (the first word of the title), so
   that the names are still unique, and the tags (tl and tr included)
   are those of the part.
*/

#ifndef TRIMMER_H
#define TRIMMER_H

#include <stdbool.h>
#include <stddef.h>
#include "amplicon.h"
#include "cleavage-site.h"
#include "markers.h"
#include "counter.h"
#include "fasta.h"
#include "fastq-reader.h"
#include "read-cover.h"
#include "reference.h"
#include "scheme.h"
#include "strbuf.h"

/* the settings of a trimmer */
typedef struct {
  const scheme_t *scheme;
  double max_errors;          /* see amplicon_analyzer_new() */
  trim_policy_t policy;
  long window;
  /* references (NULL: none), their index (built once, shared by the
     trimmers of all threads) and the settings of the classification (see
     reference.h) */
  const fasta_t *references;
  const reference_index_t *index;
  bool exhaustive;            /* see reference_set_set_exhaustive() */
  bool seeds_only;            /* see reference_set_set_seeds_only() */
  int segment;                /* see reference_set_restrict_segment() */
  size_t probe_len;
  double probe_error_rate;
  long terminal_window;
  const char *adapters;       /* the adapter kits (adapters.h; NULL or
                                 "none": none) */
  bool illumina;              /* --platform illumina: the reads are
                                 fragments of the amplicons */
  /* the anchors of the cleavage site in the HA references (NULL: not
     read; built once, shared by the trimmers of all threads) */
  const cleavage_anchors_t *cleavage;
  /* the anchors of the markers in the references (NULL: not read) */
  const marker_anchors_t *markers;
  size_t min_length;          /* shorter reads are not written */
  bool tags;                  /* append the tags to the titles */
  bool split;                 /* split the reads with several amplicons */
  bool by_category;           /* an output for each category of reads */
  bool cut_deletions;         /* cut the reads with an internal deletion
                                 at their junction */
} trimmer_config_t;

/* the greatest number of outputs of a trimmer */
#define TRIMMER_MAX_OUTPUTS READ_N_CATEGORIES

/* The number of outputs of a trimmer for <config>: one, or one for each
   category of reads (in the order of READ_CATEGORY_*) */
static inline size_t trimmer_n_outputs(const trimmer_config_t *config) {
  return config->by_category ? READ_N_CATEGORIES : 1;
}

typedef struct trimmer_t trimmer_t;

/* Create a trimmer; <config> and what it points to must exist as long as
   the trimmer is used */
trimmer_t* trimmer_new(const trimmer_config_t *config);

/*
   Trim the records <from> to <to> - 1 of <batch>; the records kept are
   appended to <out> (FASTQ; an array of trimmer_n_outputs() buffers, by
   category if more than one), the reads are counted in <counts> (created
   for the scheme and the references of the configuration)
*/
void trimmer_run(trimmer_t *self, const fastq_batch_t *batch, size_t from,
    size_t to, strbuf_t *out, counter_t *counts);

/* Classify the reads only against the references <allowed> (NULL: all);
   see amplicon_analyzer_restrict_references */
void trimmer_restrict_references(trimmer_t *self, const bool *allowed);

/* Search only the primers of the types <types>, see
   amplicon_analyzer_restrict_virus */
void trimmer_restrict_virus(trimmer_t *self, virus_set_t types);

/* Delete (NULL is accepted) */
void trimmer_delete(trimmer_t *self);

/* Unit tests */
void trimmer_test(void);

#endif
