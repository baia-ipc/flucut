/* (c) 2026 Giorgio Gonnella */

/*
   The lengths of the complete amplicons (reads with the primers of both
   ends): the length of the template between the two primers, counted in
   bins of LENGTHS_BIN bases, in groups (e.g. one for each reference, see
   trimmer.h).

   This is a check independent of the classification of each read end
   against the references: the complete amplicons of a segment have the
   length of the segment, so internal priming (or a deletion, e.g. in
   the defective interfering RNAs of the polymerase segments) shows as a
   group of shorter amplicons, at the length from the internal site to
   the other end. The check does not decide anything about a read; it is
   reported.
*/

#ifndef LENGTHS_H
#define LENGTHS_H

#include <stdio.h>
#include <stddef.h>

/* width of the bins, in bases */
#define LENGTHS_BIN 10
/* the last bin collects the lengths from LENGTHS_MAX on */
#define LENGTHS_MAX 4000
#define LENGTHS_N_BINS (LENGTHS_MAX / LENGTHS_BIN + 1)
/* an amplicon is short if its length is below this fraction of the
   expected length (e.g. that of the reference of its group) */
#define LENGTHS_SHORT_FRACTION 0.9

typedef struct lengths_t lengths_t;

/* the main length of a group and its short amplicons */
typedef struct {
  size_t total;          /* amplicons of the group */
  size_t main_bin;       /* the most frequent bin */
  size_t main_count;
  size_t short_count;    /* amplicons shorter than the fraction above of
                            the expected length */
  size_t short_bin;      /* the most frequent bin among them (if any) */
  size_t short_bin_count;
} lengths_summary_t;

/* Create the histograms of <n_groups> groups, all empty */
lengths_t* lengths_new(size_t n_groups);

/* The number of groups */
size_t lengths_n_groups(const lengths_t *self);

/* Count an amplicon of length <len> in the group <group> */
void lengths_add(lengths_t *self, size_t group, size_t len);

/* Add the counts of <other> (with the same number of groups) */
void lengths_add_all(lengths_t *self, const lengths_t *other);

/* Empty all groups */
void lengths_clear(lengths_t *self);

/* The count of bin <bin> of group <group> (the bin of a length l is
   l / LENGTHS_BIN, at most LENGTHS_N_BINS - 1) */
size_t lengths_count(const lengths_t *self, size_t group, size_t bin);

/*
   The main length of group <group>, and its short amplicons, given the
   expected length <expected> (e.g. the length of its reference: the
   main length itself is not used, as the internally primed amplicons
   can be the majority)
*/
void lengths_summarize(const lengths_t *self, size_t group,
    size_t expected, lengths_summary_t *summary);

/* a peak is a bin with at least this fraction of the amplicons of its
   group, and no larger bin within LENGTHS_PEAK_RADIUS bins */
#define LENGTHS_PEAK_FRACTION 0.02
#define LENGTHS_PEAK_RADIUS 5

/*
   The peaks of group <group> (see above), by decreasing count: at most
   <max> bins, stored in <bins>; returns their number
*/
size_t lengths_peaks(const lengths_t *self, size_t group, size_t *bins,
    size_t max);

/*
   Write the non-empty bins of group <group> to <file>, one line each:
   <prefix>, then tab-separated the start of the bin and the count (the
   prefix contains the columns before them, ending with a tab)
*/
void lengths_write(const lengths_t *self, size_t group, const char *prefix,
    FILE *file);

/* Delete (NULL is accepted) */
void lengths_delete(lengths_t *self);

/* Unit tests */
void lengths_test(void);

#endif
