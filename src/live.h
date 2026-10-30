/* (c) 2026 Giorgio Gonnella */

/*
   The live mode (--live, from the real-time mode of Lumipore): the files
   of the input directories are processed as they are written, e.g. by
   MinKNOW during a run. The directories are scanned regularly
   (--scan-interval), and the new files which are ready (not modified for
   --min-age, and complete) are trimmed; the output of each group is
   appended to, and the tables and reports are updated, after each file
   (the HTML report reloads itself while the run is in progress).

   The list of the processed files and the counts of each group are saved
   in the output directory (<group>LIVE_FILES_SFX, <group>LIVE_COUNTS_SFX),
   so that an interrupted run can be continued (--resume): the files
   already processed are skipped. Without --resume, the outputs of a
   group are replaced.

   The mode ends at the end of the run (a file LIVE_END_FILES in an input
   directory or in its parent, as MinKNOW writes it when the run ends:
   the files present are then processed, without waiting for their
   minimum age), after --timeout without new files, or when interrupted
   (Ctrl-C, SIGTERM; a second Ctrl-C kills the program); the final tables
   and reports are then written.

   -e auto is decided at the first file ready (also if already processed,
   so that a resumed run uses the same setting).
*/

#ifndef LIVE_H
#define LIVE_H

#include "auto-settings.h"
#include "fasta.h"
#include "logger.h"
#include "options.h"
#include "run-report.h"
#include "scheme.h"
#include "trimmer.h"

/* the end of the run (see above) */
#define LIVE_END_FILES "final_summary_*.txt"

/* the state files of a group (see above) */
#define LIVE_FILES_SFX ".flucut-files"
#define LIVE_COUNTS_SFX ".flucut-counts"

/*
   Run the live mode, with the options <opts>, the trimming settings
   <config> (its max_errors is set by <autos>, if pending), the <scheme>
   and the <references> (NULL: none), writing the reports with <report>
*/
void live_run(const options_t *opts, trimmer_config_t *config,
    auto_settings_t *autos, const scheme_t *scheme,
    const fasta_t *references, run_report_t *report, Logger *logger);

#endif
