/* (c) 2026 Giorgio Gonnella */

/*
   The lengths and the qualities of the reads, before and after the
   trimming: of the reads of the input, and of the reads written.

   The lengths are counted in bins of READ_STATS_LEN_BIN bases, with the
   bases of each bin, so that the median and the N50 (the length such
   that half of the bases are in reads at least that long) are those of
   the reads but for the width of a bin. The quality of a read is the
   Phred score of the mean error probability of its bases (as in MinKNOW
   and NanoPlot: not the mean of the scores, which would hide the bad
   bases), counted in bins of one unit. The memory does not depend on
   the number of reads.
*/

#ifndef READ_STATS_H
#define READ_STATS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

/* width of the bins of the lengths; the last bin collects the lengths
   from READ_STATS_LEN_MAX on */
#define READ_STATS_LEN_BIN 20
#define READ_STATS_LEN_MAX 20000
#define READ_STATS_LEN_BINS (READ_STATS_LEN_MAX / READ_STATS_LEN_BIN + 1)
/* the qualities from this value on are counted together */
#define READ_STATS_QUAL_MAX 60
#define READ_STATS_QUAL_BINS (READ_STATS_QUAL_MAX + 1)

/* Maximum number of independently counted read sets in one object. */
#define READ_STATS_MAX_SETS 8

typedef struct {
  size_t reads, bases;
  size_t min_len, max_len;         /* (0 without reads) */
  size_t len_reads[READ_STATS_LEN_BINS];
  size_t len_bases[READ_STATS_LEN_BINS];
  size_t with_quality;             /* reads with qualities (not empty) */
  size_t quality[READ_STATS_QUAL_BINS];
} read_stats_set_t;

typedef struct {
  size_t n_sets;
  const char *const *set_names; /* caller-owned, must outlive this object */
  read_stats_set_t sets[READ_STATS_MAX_SETS];
} read_stats_t;

/* what the reports say of a set of reads */
typedef struct {
  size_t reads, bases;
  size_t min_len, max_len, mean_len, median_len, n50;
  bool has_quality;                /* false: the fields below are 0 */
  double median_quality;
  double q10, q20;                 /* shares of the reads with a quality
                                      of at least 10, of at least 20 */
} read_stats_summary_t;

/* Initialize with <n_sets> names; names and their strings stay caller-owned. */
void read_stats_init(read_stats_t *self, size_t n_sets,
    const char *const *set_names);
/* Empty the counts while keeping the set configuration. */
void read_stats_clear(read_stats_t *self);

/* Count a read of <len> bases in the set <set>; <qual>: its qualities
   (<len> characters, Phred+33), or NULL if it has none */
void read_stats_add(read_stats_t *self, size_t set, size_t len,
    const char *qual);

/* Add the counts of <other> */
void read_stats_merge(read_stats_t *self, const read_stats_t *other);

/* The summary of the set <set> */
void read_stats_summarize(const read_stats_t *self, size_t set,
    read_stats_summary_t *out);

/* The configured name of a set. */
const char* read_stats_set_name(const read_stats_t *self, size_t set);

/*
   Save the counts to <file>, as lines of text starting with "rlen",
   "rqual" or "rset"; and load such a line, added to the counts: false
   if <line> is not a valid one
*/
void read_stats_save(const read_stats_t *self, FILE *file);
bool read_stats_load_line(read_stats_t *self, const char *line);

/* Unit tests */
void read_stats_test(void);

#endif
