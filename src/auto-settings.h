/* (c) 2026 Giorgio Gonnella */

/*
   The automatic settings (from Lumipore): -e auto, the errors allowed in
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

/* Free the memory */
void auto_settings_free(auto_settings_t *self);

/* Unit tests */
void auto_settings_test(void);

#endif
