/* (c) 2026 Giorgio Gonnella */

/*
   Identification of the cell barcodes of a single-cell library: the
   barcodes of a whitelist (e.g. all barcodes of a platform) found exactly
   in the reads are counted, and the barcodes of the cells are selected
   from the counts, given the expected number of cells N, as in Cell
   Ranger (version 2): the barcodes with a count of at least 1/10 of the
   99th percentile of the counts of the N most frequent barcodes; as in
   Cell Ranger, the threshold is rounded to the nearest integer (halves
   upwards; e.g. 14 for a percentile of 135..144), and is at least 1.
*/

#ifndef CELL_CALLER_H
#define CELL_CALLER_H

#include <stddef.h>
#include <stdint.h>
#include "logger.h"
#include "whitelist.h"

typedef struct cell_caller_t cell_caller_t;

cell_caller_t* cell_caller_new(void);

/* Count an occurrence of an encoded barcode */
void cell_caller_add(cell_caller_t *self, uint64_t barcode);

/* Add the counts of <other> and reset them */
void cell_caller_merge(cell_caller_t *self, cell_caller_t *other);

/*
   Select the cell barcodes, given the expected number of cells; returns
   them as a whitelist (to be deleted by the caller), with barcodes of
   length <length>; <known>: if not NULL, cells known from another source
   (e.g. the cells called by Cell Ranger in an Illumina run of the same
   library), which are added to the selected ones (e.g. for targeted
   libraries with few reads per cell, where only the cells with most reads
   pass the threshold); if <path> is not NULL, the cells are written to it,
   with their counts (the selected ones, then the known ones not selected)
*/
whitelist_t* cell_caller_run(const cell_caller_t *self,
    size_t expected_cells, size_t length, const whitelist_t *known,
    const char *path, Logger *logger);

/*
   Counts of the barcodes, in decreasing order (allocated, to be freed by
   the caller); returns their number
*/
size_t cell_caller_sorted_counts(const cell_caller_t *self,
    uint64_t **counts);

/* Delete a cell caller (NULL is accepted) */
void cell_caller_delete(cell_caller_t *self);

/* Unit tests */
void cell_caller_test(const char *tmpdir);

#endif
