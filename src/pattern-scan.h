/* (c) 2026 Giorgio Gonnella */

/*
   All the places of a sequence where a short pattern ends, within a
   number of errors, in one pass over the sequence.

   The search of a primer with edlib (primer-search.h) gives the best
   match of a pattern in a sequence; when all the matches are needed
   (the termini inside a read, see chimera.h; the primers inside a read,
   see amplicon.h), the sequence would have to be searched again on both
   sides of each match. Here the sequence is read once: the edit distance
   of the pattern to the bases which end at each position is kept up to
   date by the bit-parallel algorithm of Myers (1999, J ACM 46:395-415,
   https://doi.org/10.1145/316542.316550), one machine word for a
   pattern of up to PATTERN_SCAN_WORD bases; a longer pattern is handled
   by the plain dynamic programming (one column of its matrix).

   The characters are compared as in the search of a primer (see
   iupac.h).

   The scan tells where a match ends, and with how many errors
   (pattern_scan_run(), pattern_scan_best()); for the matches themselves
   (pattern_scan_matches()) the pattern is then aligned backwards from
   each end to find where the match begins.

   A prepared pattern contains buffers: use one for each thread.
*/

#ifndef PATTERN_SCAN_H
#define PATTERN_SCAN_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* the longest pattern scanned with one machine word */
#define PATTERN_SCAN_WORD 64

/* a pattern prepared for the scan */
typedef struct pattern_scan_t pattern_scan_t;

/*
   Consecutive positions (0-based, <first> to <last>) of a sequence at
   which the pattern ends with at most the errors allowed: one match
   usually ends acceptably at a few positions next to each other.
   <best>: the last of them with the fewest errors, <n_errors>.
*/
typedef struct {
  long first, last, best;
  long n_errors;
} pattern_scan_end_t;

/* Prepare the pattern <pattern> (length <len>, at least 1) */
pattern_scan_t* pattern_scan_new(const char *pattern, size_t len);

/* Replace the pattern of <self> by <pattern> (length <len>, at least 1) */
void pattern_scan_set(pattern_scan_t *self, const char *pattern,
    size_t len);

/*
   Scan <seq> (length <len>) for the pattern with at most <max_errors>
   errors (fewer than the length of the pattern): the runs of positions
   at which it ends are stored in <ends> (room for <max_ends>), in the
   order of the sequence, and their number is returned (the runs after
   the first <max_ends> are not reported).
*/
size_t pattern_scan_run(pattern_scan_t *self, const char *seq, size_t len,
    long max_errors, pattern_scan_end_t *ends, size_t max_ends);

/*
   The fewest errors with which the pattern ends in <seq> (length
   <len>), if at most <max_errors> (else false is returned): <n_errors>,
   and the first and the last position at which it ends with that many,
   <first> and <last>. <backwards>: <seq> is read from its last base to
   its first, and the positions are those at which the pattern BEGINS;
   the pattern must then be given reversed (its last base first).
*/
bool pattern_scan_best(pattern_scan_t *self, const char *seq, size_t len,
    long max_errors, bool backwards, long *n_errors, long *first,
    long *last);

/* the result of pattern_scan_best(), for one of two patterns */
typedef struct {
  bool found;
  long n_errors, first, last;
} pattern_scan_best_t;

/*
   Two patterns in one pass over <seq> (length <len>): <backwards> as by
   pattern_scan_best() reading <seq> from its last base to its first
   (into <b>), <forwards> reading it from its first to its last (into
   <f>), both with at most <max_errors> errors. The results are those of
   the two calls of pattern_scan_best(); for two patterns of at most
   PATTERN_SCAN_WORD bases it takes less time than the two calls.
*/
void pattern_scan_best_both(pattern_scan_t *backwards,
    pattern_scan_t *forwards, const char *seq, size_t len, long max_errors,
    pattern_scan_best_t *b, pattern_scan_best_t *f);

/* the most matches of a pattern reported in one sequence */
#define PATTERN_SCAN_MAX_MATCHES 32

/* a match of a pattern (0-based, inclusive) */
typedef struct {
  long start, end;
  long n_errors;
} pattern_match_t;

/*
   The matches of the pattern in <seq> (length <len>) with at most
   <max_errors> errors, in the order of the sequence, into <matches>
   (room for PATTERN_SCAN_MAX_MATCHES); their number is returned. A
   match ends at the last of the consecutive positions with the fewest
   errors, and begins as early as that number of errors allows (the
   longest match). Matches of the same pattern can overlap (e.g. in a
   repeat).
*/
size_t pattern_scan_matches(pattern_scan_t *self, const char *seq,
    size_t len, long max_errors, pattern_match_t *matches);

/* Delete a prepared pattern (NULL is accepted) */
void pattern_scan_delete(pattern_scan_t *self);

#endif
