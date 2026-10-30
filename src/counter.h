/* (c) 2026 Giorgio Gonnella */

/*
   The counts of the reads of a group of input files (or of all reads),
   for the summary table and the reports:
   - the reads at each step (counter_reads_t);
   - the reads with a primer found at each end, by primer of the scheme
     (variant);
   - the internal priming sites: with references, the read ends
     classified as internal, by reference, end and position (that of the
     first base after the primer, see reference_hit_t);
   - the lengths of the complete amplicons (lengths.h).

   A counter belongs to one thread; the counters of the threads are
   summed (counter_merge()) in the order of the input. The memory does
   not depend on the number of reads: the sites are counted in an array
   of the positions of each reference, allocated when the reference has
   its first internal end (at most twice the total length of the
   references).
*/

#ifndef COUNTER_H
#define COUNTER_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include "amplicon.h"
#include "fasta.h"
#include "lengths.h"
#include "scheme.h"

/* the ends of an amplicon */
enum { COUNTER_LEFT, COUNTER_RIGHT, COUNTER_N_ENDS };

/* the counts of the reads at each step */
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
  /* reads with an adapter removed at their start, at their end (see
     amplicon.h) */
  size_t adapter_start, adapter_end;
  /* ends found, by end (COUNTER_LEFT, COUNTER_RIGHT) and origin
     (end_origin_t); all unknown without references */
  size_t origin[COUNTER_N_ENDS][3];
} counter_reads_t;

typedef struct counter_t counter_t;

/*
   Create an empty counter for the primers of <scheme> and the
   <references> (NULL: none); both must exist as long as the counter is
   used
*/
counter_t* counter_new(const scheme_t *scheme, const fasta_t *references);

/* Empty the counter */
void counter_clear(counter_t *self);

/* Add the counts of <other> (created with the same scheme and
   references) */
void counter_merge(counter_t *self, const counter_t *other);

/*
   Count the analysis <a> of a read of length <len> (the reads written or
   too short are counted by counter_add_written() or
   counter_add_too_short())
*/
void counter_add_read(counter_t *self, const amplicon_t *a, size_t len);

/* Count a read written, with <len> bases */
void counter_add_written(counter_t *self, size_t len);

/* Count a read not written, as too short */
void counter_add_too_short(counter_t *self);

/* The counts of the reads at each step */
const counter_reads_t* counter_reads(const counter_t *self);

/* The number of primers of an end (COUNTER_LEFT or COUNTER_RIGHT) of the
   scheme */
size_t counter_n_variants(const counter_t *self, size_t end);

/* The reads in which the primer <variant> of the scheme was found at the
   end <end> */
size_t counter_primer(const counter_t *self, size_t end, size_t variant);

/*
   The internal ends of the reference <ref>, at the end <end>: an array
   of the length of the reference, by 0-based position on its plus
   strand (of the first base after the primer, see reference_hit_t);
   NULL if there are none (or no references)
*/
const size_t* counter_sites(const counter_t *self, size_t ref, size_t end);

/*
   The lengths of the complete amplicons; the groups: without references,
   one; with references, one for each reference (the reference to which
   the probe of the left end, or else of the right end, mapped) and a
   last one for the amplicons not located
*/
const lengths_t* counter_lengths(const counter_t *self);

/* The number of groups of the lengths of a counter for <references>
   (NULL: none), see counter_lengths() */
size_t counter_n_length_groups(const fasta_t *references);

/*
   Save the counts to <file> (text; for the resumption of a run in
   progress, see live.h), and load them, added to those of <self> (created
   with the same scheme and references); false if <file> is not a valid
   file of counts (then the counts of <self> are undefined)
*/
void counter_save(const counter_t *self, FILE *file);
bool counter_load(counter_t *self, FILE *file);

/* Delete (NULL is accepted) */
void counter_delete(counter_t *self);

/* Unit tests */
void counter_test(void);

#endif
