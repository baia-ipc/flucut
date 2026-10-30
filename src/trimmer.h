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

   With references, for each end found:

     ol:Z:, or:Z:   origin of the left/right end: terminal, internal or
                    unknown (see amplicon.h)
     ml:Z:, mr:Z:   where its probe mapped: <reference>:<pos>:<strand>
                    (1-based position on the plus strand)

   Reads shorter than the minimum length after the trimming are not
   written (primer dimers, with nothing left, always).
*/

#ifndef TRIMMER_H
#define TRIMMER_H

#include <stdbool.h>
#include <stddef.h>
#include "amplicon.h"
#include "fasta.h"
#include "fastq-reader.h"
#include "lengths.h"
#include "reference.h"
#include "scheme.h"
#include "strbuf.h"

/* the counts of the reads */
typedef struct {
  size_t reads_in, bases_in;
  size_t reads_out, bases_out;
  /* reads not written, as shorter than the minimum length after the
     trimming (dimers included) */
  size_t too_short;
  /* reads with both ends found and nothing between them */
  size_t dimers;
  /* reads by number of ends found (0, 1, 2) */
  size_t n_ends[3];
  /* reads whose layout was found on the minus strand */
  size_t minus;
  /* reads with primer matches inside them */
  size_t with_internal;
  /* ends found, by end (0: left, 1: right) and origin (end_origin_t);
     all unknown without references */
  size_t origin[2][3];
} trimmer_stats_t;

/* Add the counts <other> to <self> */
void trimmer_stats_add(trimmer_stats_t *self, const trimmer_stats_t *other);

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
  size_t probe_len;
  double probe_error_rate;
  long terminal_window;
  size_t min_length;          /* shorter reads are not written */
  bool tags;                  /* append the tags to the titles */
} trimmer_config_t;

typedef struct trimmer_t trimmer_t;

/* Create a trimmer; <config> and what it points to must exist as long as
   the trimmer is used */
trimmer_t* trimmer_new(const trimmer_config_t *config);

/*
   Trim the records <from> to <to> - 1 of <batch>; the records kept are
   appended to <out> (FASTQ), the counts added to <stats>, the lengths of
   the complete amplicons (the template between the two primers, primer
   dimers excluded) to <lengths> (see lengths.h; its groups: see
   trimmer_n_length_groups())
*/
void trimmer_run(trimmer_t *self, const fastq_batch_t *batch, size_t from,
    size_t to, strbuf_t *out, trimmer_stats_t *stats, lengths_t *lengths);

/*
   The groups of the lengths: without references, one; with references,
   one for each reference (the reference to which the probe of the left
   end, or else of the right end, mapped) and a last one for the
   amplicons not located
*/
size_t trimmer_n_length_groups(const trimmer_config_t *config);

/* Delete (NULL is accepted) */
void trimmer_delete(trimmer_t *self);

/* Unit tests */
void trimmer_test(void);

#endif
