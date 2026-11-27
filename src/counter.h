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
   - the lengths of the complete amplicons (lengths.h);
   - the lengths and qualities of the reads, of the input and written
     (read-stats.h), and what their titles say about the sequencing run
     (read-metadata.h).

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
#include "cleavage-site.h"
#include "fasta.h"
#include "lengths.h"
#include "read-cover.h"
#include "read-metadata.h"
#include "read-stats.h"
#include "scheme.h"

/* the ends of an amplicon */
enum { COUNTER_LEFT, COUNTER_RIGHT, COUNTER_N_ENDS };

/* the counts of the reads at each step */
typedef struct {
  size_t reads_in, bases_in;
  size_t reads_out, bases_out;
  /* the reads written, by category (READ_CATEGORY_*, see read-cover.h);
     without references all are READ_CATEGORY_OTHER */
  size_t written[READ_N_CATEGORIES];
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
  /* reads split into several amplicons (chimera.h); each part is then
     counted as a read in all the counts but <reads_in> and <bases_in>,
     which are those of the input */
  size_t split;
  /* reads with an adapter removed at their start, at their end (see
     amplicon.h) */
  size_t adapter_start, adapter_end;
  /* ends found, by end (COUNTER_LEFT, COUNTER_RIGHT) and origin
     (end_origin_t); all unknown without references */
  size_t origin[COUNTER_N_ENDS][3];
} counter_reads_t;

/* The reads analyzed: those of the input, each part of a split read
   counted as one (the reads to which the counts other than <reads_in>
   and <bases_in> refer) */
static inline size_t counter_reads_analyzed(const counter_reads_t *s) {
  return s->n_ends[0] + s->n_ends[1] + s->n_ends[2];
}

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

/* The same, for a read already located on the references of the counter
   (<cover>, see read_cover_locate()) */
void counter_add_located_read(counter_t *self, const amplicon_t *a,
    size_t len, const read_cover_t *cover);

/*
   Count a read of length <len> split into <n_parts> parts (chimera.h),
   of <parts_len> bases together, each of which was counted by
   counter_add_read(): the reads and the bases of the input are counted
   once
*/
void counter_add_split(counter_t *self, size_t n_parts, size_t len,
    size_t parts_len);

/*
   Count a read of the input for the statistics of the reads: its <len>
   bases, with the qualities <qual> (NULL: none); <title>: its title, of
   length <title_len>, read for the metadata of the run (NULL: not read;
   the titles of some reads are enough, as those of a run say the same)
*/
void counter_add_input(counter_t *self, size_t len, const char *qual,
    const char *title, size_t title_len);

/* Count a read written, with <len> bases, of the category <category>
   (READ_CATEGORY_*), with the qualities <qual> (NULL: none) */
void counter_add_written(counter_t *self, size_t len, size_t category,
    const char *qual);

/*
   Count a read the cleavage site of the haemagglutinin was looked for
   in, and the motif <nt> read from one (cleavage-site.h); which reads
   are read is decided by the trimmer
*/
void counter_set_cleavage_reference(counter_t *self, const char *aa);
void counter_add_cleavage_searched(counter_t *self);
void counter_add_cleavage(counter_t *self, const char *nt);

/* The motifs of the cleavage site read from the reads */
const cleavage_counts_t* counter_cleavage(const counter_t *self);

/* The lengths and qualities of the reads, and the metadata of the run */
const read_stats_t* counter_read_stats(const counter_t *self);
const read_metadata_t* counter_metadata(const counter_t *self);

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
   The reconstructed coverage of the reference <ref> by the class <cat>
   (COVER_*) of the located reads, from their spans (see read-cover.h):
   <out> is filled with the number of reads of that class covering each
   of the <seq_len> positions of the reference, and the greatest value
   is returned (0 and <out> zeroed if none). <out> must have room for
   the reference length.
*/
size_t counter_coverage(const counter_t *self, size_t ref, size_t cat,
    size_t *out);

/* Whether any located read was counted for the coverage of <ref> */
bool counter_has_coverage(const counter_t *self, size_t ref);

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
