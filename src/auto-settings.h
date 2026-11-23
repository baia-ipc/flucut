/* (c) 2026 Giorgio Gonnella */

/*
   The automatic settings: -e auto, the errors allowed in
   a primer, decided from the error rate of the primer alignments in the
   first AUTO_MAX_ERRORS_READS reads of the first input file, before the
   trimming (in the live mode, at the first file ready). The decision and
   the measurement are kept, for the reports (HTML report, JSON summary).

   The error rate is that of the alignments of the primers found at the
   ends of the reads with AUTO_MAX_ERRORS_LOW (errors over aligned read
   bases). The errors allowed are then AUTO_MAX_ERRORS_LOW, or
   AUTO_MAX_ERRORS_HIGH from a rate of AUTO_MAX_ERRORS_HIGH_RATE: on reads
   simulated with mean error rates of 1, 2, 4, 8 and 12%
   (benchmark/RESULTS.md), whose measured rates were 0.9, 1.7, 3.0, 5.0
   and 6.3%, 0.2 (2 errors in Uni12) left few primer bases up to 4%,
   while at 8 and 12% 0.3 (3 errors) left 4 times fewer primer bases and
   2 to 4 times fewer bases of the internal priming site than 0.2, at the
   cost of 3 times more template bases removed.
*/

#ifndef AUTO_SETTINGS_H
#define AUTO_SETTINGS_H

#include <stddef.h>
#include "logger.h"
#include "options.h"
#include "report.h"
#include "trimmer.h"

#define AUTO_MAX_ERRORS_READS 1000
/* the choice of the reference of each segment for a group (see
   auto_settings_choose_references): the reads used at first, and in all
   if a segment has too few votes in them; the probes taken from each
   read; the fewest votes of the reference chosen; its share of the
   votes of its segment from which its lead is clear; the share of the
   votes of the leading reference from which another one is a close
   contender */
#define AUTO_REFERENCE_READS 2000
#define AUTO_REFERENCE_MORE_READS 20000
#define AUTO_REFERENCE_PROBES 8
#define AUTO_REFERENCE_MIN_VOTES 20
#define AUTO_REFERENCE_MIN_SHARE 0.6
/* (on simulated reads of 17 genomes outside the built-in set, the three
   values above change the share of read ends classified correctly by
   less than a point: above a close share of 0.25 the reference of the
   common subtype and the leading one are equally good, below it the
   first is worse; benchmark/RESULTS.md) */
#define AUTO_REFERENCE_CLOSE_SHARE 0.5
/* (the help text of -e, in options.c, states these three values) */
#define AUTO_MAX_ERRORS_LOW 0.2
#define AUTO_MAX_ERRORS_HIGH 0.3
#define AUTO_MAX_ERRORS_HIGH_RATE 0.04

typedef struct {
  /* -e: the value, and if automatic the measurement (see report.h) */
  report_setting_t max_errors;
  char *file;   /* the file of the measurement (owned) */
} auto_settings_t;

/* Initialize the settings from the options <opts> */
void auto_settings_init(auto_settings_t *self, const options_t *opts);

/* Whether a setting is still to be decided */
bool auto_settings_pending(const auto_settings_t *self);

/* The errors allowed in a primer for the error rate <error_rate> */
double auto_settings_max_errors(double error_rate);

/*
   The error rate of the primer alignments in the first <max_reads> reads
   of <file>, with the settings of <config> but AUTO_MAX_ERRORS_LOW errors
   allowed (see above); the number of reads used is stored in <n_reads>
*/
double auto_settings_error_rate(const trimmer_config_t *config,
    const char *file, size_t max_reads, size_t *n_reads, Logger *logger);

/*
   Decide the settings still to be decided from the first reads of <file>,
   and set them in <config> (the decisions are reported to <logger>)
*/
void auto_settings_decide(auto_settings_t *self, trimmer_config_t *config,
    const char *file, Logger *logger);

/*
   Choose, for each segment, the reference closest to the reads of
   <file>, among the references of that segment (e.g. those of several
   subtypes). Each group (a barcode, a different sample) is then
   classified against the chosen references only, so that the other
   references of a segment do not match it spuriously. Each segment is
   chosen by its own reads, as the segments of a virus can be closest
   to references of different subtypes (a reassortant); a segment whose
   reads do not single out a reference follows the clear ones.

   Probes are taken along the first AUTO_REFERENCE_READS reads
   (AUTO_REFERENCE_PROBES from each, evenly spaced) and located by their
   k-mers; a probe votes for the reference it is located in, unless
   another reference matches it as well (it then tells nothing).

   A segment is clear if its reference with most votes has at least
   AUTO_REFERENCE_MIN_VOTES votes and at least AUTO_REFERENCE_MIN_SHARE
   of the votes of the segment: that reference is chosen. The other
   segments follow the clear ones, where they can: the subtype chosen
   in most clear segments (in case of a tie, the one with most votes in
   them) is the common subtype, and its reference is chosen
   - in a segment with too few votes;
   - in a segment without a clear lead, if it is a close contender (at
     least AUTO_REFERENCE_CLOSE_SHARE of the votes of the leading
     reference).
   Where this does not apply (no clear segment, references without a
   subtype, the reference of the common subtype is not close), a
   segment without a clear lead takes its reference with most votes
   (the sample is about as close to other references, and more reads
   would not change that), and for a segment with too few votes the
   reads up to AUTO_REFERENCE_MORE_READS are used as well; if it still
   has too few, it keeps all its references. The choices are reported.

   Returns a new array with a value for each reference (whether it is
   considered, see reference_set_restrict), to be freed by the caller,
   or NULL if all references are considered (no segment has several
   references, or none was decided).
*/
bool* auto_settings_choose_references(const trimmer_config_t *config,
    const char *file, Logger *logger);

/*
   The decision of auto_settings_choose_references(), from the votes:
   <segment>, <subtype> and <votes> have a value for each of the <n>
   references (<subtype> is an id, -1 if not known; NULL if none is
   known); <allowed> is set for each. Returns the number of segments
   with several references and too few votes for a choice.
*/
size_t auto_settings_decide_references(const int *segment,
    const int *subtype, const size_t *votes, size_t n, bool *allowed);

/*
   The common subtype of the clear segments (see
   auto_settings_choose_references()), -1 if there is none
*/
int auto_settings_common_subtype(const int *segment, const int *subtype,
    const size_t *votes, size_t n);

/* Whether a reference with <best_votes> of the <segment_votes> votes of
   its segment has a clear lead */
bool auto_settings_clear_lead(size_t best_votes, size_t segment_votes);

/* Free the memory */
void auto_settings_free(auto_settings_t *self);

/* Unit tests */
void auto_settings_test(void);

#endif
