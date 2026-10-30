/* (c) 2017-2021 Giorgio Gonnella */
/* Merging of sorted files */

#ifndef FILE_MERGER_H
#define FILE_MERGER_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include "logger.h"

/*
  Function reading an element: returns true if an element was obtained from
  file, false if eof was reached and exits the program if an error occurred;
  furthermore, it stores the element in a char* buffer, which is pointed by
  elem, allocated if NULL and reallocated if not large enough; its size is
  stored in the size_t pointed by elemsize (0 at the first call for a
  file, when the buffer is NULL).
*/
typedef bool (*merge_getelem_f)(char **elem, size_t *elemsize, FILE *file,
    Logger *logger);

/* Function writing the element elem to file */
typedef void (*merge_writeelem_f)(FILE *file, char *elem, Logger *logger);

/*
  Compare function (as for qsort); the arguments are pointers to the
  element buffers, i.e. char ** (cast from const void *), not the elements
*/
typedef int (*merge_cmpelems_f)(const void *elem1, const void *elem2);

/*
  Merge n text files element-wise using a compare function.
  The input files are named <basename><sfx>.N (e.g. with sfx ".sorted":
  big.txt.sorted.0), where N is a number between 0 and nfiles-1, as
  written by sort_file_parts() of file-sorter.h; the input files are
  deleted. At most 512 files are open at the same time: more files are
  merged in batches, into further temporary files <basename><sfx>.N
  (N from nfiles on), which are then merged in the same way.

  Parameters:
  - nfiles: number of files to merge
  - basename: basename of the input/output files
  - sfx: suffix of the output and temporary files (appended to basename,
         thus usually starting with a dot)
  - use_stdout: if true, output is to the standard output;
                otherwise to <basename><sfx>
  - getelem, writeelem, cmpelems: see the types above
*/
void merge_n_files(size_t nfiles, const char *basename, const char *sfx,
    bool use_stdout, merge_getelem_f getelem, merge_writeelem_f writeelem,
    merge_cmpelems_f cmpelems, Logger *logger);

#endif
