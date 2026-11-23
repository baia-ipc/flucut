/* (c) 2026 Giorgio Gonnella */

/*
   The state of a group (a sample), in the subdirectory GROUP_STATE_DIR
   of the output directory, so that a run can be continued (--resume):
   the counts of the reads processed (<group>GROUP_COUNTS_SFX) and the
   list of the input files processed (<group>GROUP_FILES_SFX).

   Each file is written to a temporary file and then renamed, so that an
   interruption while writing it does not leave it incomplete. The list
   is written after the counts: a file is only listed once its counts are
   saved, so after an interruption between the two the group (or, in the
   live mode, the file) is processed again, instead of being skipped with
   its counts missing.

   In the live mode (live.h) the list grows file by file. Otherwise it is
   written once a group is complete, and kept after the end of the run:
   --resume on a complete run then only writes the reports again.
*/

#ifndef GROUP_STATE_H
#define GROUP_STATE_H

#include <stdbool.h>
#include <stddef.h>
#include "counter.h"
#include "logger.h"

#define GROUP_STATE_DIR ".flucut-state"
#define GROUP_FILES_SFX ".flucut-files"
#define GROUP_COUNTS_SFX ".flucut-counts"

/* Save the counts of the group <group> of <outdir>, then the list of its
   <n_files> processed <files> (GROUP_STATE_DIR is created if needed) */
void group_state_save(const char *outdir, const char *group,
    const counter_t *counts, char *const *files, size_t n_files,
    Logger *logger);

/* The list of the processed files of a group (to be freed with
   group_state_free_files()); NULL, with <n_files> 0, if there is none */
char** group_state_load_files(const char *outdir, const char *group,
    size_t *n_files);
void group_state_free_files(char **files, size_t n_files);

/* Add the saved counts of a group to <counts>; false if there are none,
   or they cannot be read */
bool group_state_load_counts(const char *outdir, const char *group,
    counter_t *counts);

/* Remove the state files of a group, if any */
void group_state_remove(const char *outdir, const char *group,
    Logger *logger);

#endif
