/* (c) 2017-2021 Giorgio Gonnella */

/* Sorting of large files in parts which fit in memory (to be merged with
   merge_n_files() of file-merger.h) */

#ifndef FILE_SORTER_H
#define FILE_SORTER_H

#include <stddef.h>
#include <stdio.h>
#include "logger.h"

/*
  Read the input file in parts of up to <bufsize> bytes, split each part
  into elements (read_elems(), e.g. read_lines_in_buffer() of file.h), sort
  them (sort_elems()) and output them (output_sorted_elems()) to the file
  <outfname><outfsfx>.N, where N is the number of the part.

  Returns the number of parts.
*/
size_t sort_file_parts(
    const char *infname,  /* name (or basename) of input file */
    const char *infsfx,   /* suffix of input file (use NULL if n.a.) */
    const char *outfname, /* basename for temp and output files */
    const char *outfsfx,  /* sfx for temp (before num) and output files */
    size_t bufsize,       /* size of memory buffer */
    size_t elems_bufsize_init, /* initial size of the array of the
                                  element pointers, in bytes, at least 1
                                  (it is enlarged by read_elems() as
                                  needed) */
    size_t read_elems(FILE *file, char *buffer, size_t bufsize,
                      char ***elems, size_t *elemsbufsize, Logger *logger),
    void sort_elems(char **elems, size_t nelems),
    void output_sorted_elems(FILE *file, char **elems, size_t nelems,
                             Logger *logger),
    Logger *logger);

#endif
