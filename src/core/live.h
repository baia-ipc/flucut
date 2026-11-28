/* (c) 2026 Giorgio Gonnella */

/*
   The live mode (--live): the files
   of the input directories are processed as they are written, e.g. by
   MinKNOW during a run. The directories are scanned regularly
   (--scan-interval), and the new files which are ready (not modified for
   --min-age, and complete) are trimmed; the output of each group is
   appended to, and the tables and reports are updated, after each file
   (the HTML report reloads itself while the run is in progress).

   The list of the processed files and the counts of each group are saved
   in the output directory after each file (group-state.h), so that an
   interrupted run can be continued (--resume): the files already
   processed are skipped. Without --resume, the outputs of a
   group are replaced.

   The mode ends at the end of the run (a file LIVE_END_FILES in an input
   directory or in its parent, as MinKNOW writes it when the run ends:
   the files present are then processed, without waiting for their
   minimum age), after --timeout without new files, or when interrupted
   (Ctrl-C, SIGTERM; a second Ctrl-C kills the program); the final tables
   and reports are then written.

   -e auto is decided at the first file ready (also if already processed,
   so that a resumed run uses the same setting). With -s auto the
   scheme is chosen from that file as well, before anything is counted
   or reported with it: the mode waits for it first.
*/

#ifndef LIVE_H
#define LIVE_H

#include "auto-settings.h"
#include "fasta.h"
#include "logger.h"
#include "options.h"
#include "platform-check.h"
#include "run-report.h"
#include "scheme.h"
#include "scheme-choice.h"
#include "trimmer.h"

/* the end of the run (see above) */
#define LIVE_END_FILES "final_summary_*.txt"

/*
   Run the live mode, with the options <opts>, the trimming settings
   <config> (those of the platform set by <platform>, and its max_errors
   by <autos>, if pending), the <scheme> (chosen by <choice>, if
   pending) and the <references> (NULL: none), writing the reports with
   <report>
*/
void live_run(options_t *opts, trimmer_config_t *config,
    platform_check_t *platform, scheme_choice_t *choice,
    auto_settings_t *autos, scheme_t *scheme, const fasta_t *references,
    run_report_t *report, Logger *logger);

#endif
