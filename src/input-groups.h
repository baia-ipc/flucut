/* (c) 2026 Giorgio Gonnella */

/*
   Groups of input files (from Lumipore): each group gives one output.
   - a file of reads given directly is a group, named as the file without
     its suffix (.fastq, .fq, optionally .gz, or .bam);
   - a directory is a group of its files of reads, named as the
     directory, unless it contains subdirectories whose name starts with
     INPUT_GROUPS_DIR_PFX (e.g. the barcode directories of the fastq_pass
     directory of MinKNOW): then each of them is a group, named as the
     subdirectory, and the files beside them are not read (a warning is
     output).
   The output of a group is <outdir>/<group>.fastq.gz
   (input_group_output_path()).
   input_groups_*(): the array of groups; input_group_*(): a single group.
*/

#ifndef INPUT_GROUPS_H
#define INPUT_GROUPS_H

#include <stdbool.h>
#include <stddef.h>
#include "logger.h"

#define INPUT_GROUPS_DIR_PFX "barcode"
/* suffix of the output file of a group */
#define INPUT_GROUPS_OUTPUT_SFX ".fastq.gz"

typedef struct {
  char *name;
  char **files;    /* sorted by name */
  size_t n_files;
} input_group_t;

/*
   The groups of the <n_inputs> <inputs> (files or directories), whose
   number is stored in <n_groups>. A warning is output for a group
   without files, unless <allow_empty> (e.g. a run in progress: the files
   are still to come), and, if <warn_ignored>, for the files of reads of
   a directory beside its barcode directories (false for repeated scans,
   which would repeat it). The program is exited with an error if a file
   given directly is not a file of reads (by its suffix) or if the names
   of two groups are equal (as one output would overwrite the other).
*/
input_group_t *input_groups_new(char *const *inputs, size_t n_inputs,
    size_t *n_groups, bool allow_empty, bool warn_ignored, Logger *logger);

/* Delete the groups (NULL is accepted, if <n_groups> is 0) */
void input_groups_delete(input_group_t *groups, size_t n_groups);

/* Path of the output file of the group <group> (to be freed by the
   caller) */
char *input_group_output_path(const char *outdir, const char *group);

/* Unit tests (<tmpdir>: directory for temporary files, left empty) */
void input_groups_test(const char *tmpdir);

#endif
