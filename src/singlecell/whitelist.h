/* (c) 2026 Giorgio Gonnella */

/*
   Whitelist of barcodes (e.g. the cell barcodes of a single-cell
   platform) and correction of the barcodes extracted from the reads.

   The whitelist is read from a text file (optionally gzipped), with one
   barcode per line (only the first word of each line is used, without a
   suffix starting with "-", such as the "-1" in the barcode files of
   Cell Ranger); all barcodes must have the same length (at most 32) and
   consist of the letters ACGT.

   A barcode is searched at one end of a sequence (the anchored end, e.g.
   the end adjacent to an adaptor), allowing up to a maximum edit distance:
   the whitelist barcode must be aligned completely, while the sequence
   can be aligned partially (i.e. the other end is free, as the barcode is
   followed by other elements, e.g. an UMI). The barcode with the minimum
   distance is returned if it is unique; for each distance, all
   sequences at that distance are enumerated and looked up, thus the time
   increases quickly with the distance (a distance of 2 is the practical
   limit).
*/

#ifndef WHITELIST_H
#define WHITELIST_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "logger.h"

/* maximum length of the barcodes */
#define WHITELIST_MAX_BARCODE_LENGTH 32

/* maximum supported edit distance */
#define WHITELIST_MAX_DIST 3

typedef struct whitelist_t whitelist_t;

typedef enum {
  WHITELIST_EXACT,      /* found with distance 0 */
  WHITELIST_CORRECTED,  /* unique barcode with the minimum distance > 0 */
  WHITELIST_AMBIGUOUS,  /* several barcodes with the minimum distance */
  WHITELIST_NOT_FOUND,  /* no barcode within the maximum distance */
} whitelist_status_t;

typedef struct {
  whitelist_status_t status;
  /* exact and corrected only: */
  uint64_t barcode;     /* encoded barcode, see whitelist_decode() */
  long dist;            /* edit distance */
  size_t consumed;      /* length of the aligned part of the sequence */
} whitelist_match_t;

/*
   Read a whitelist from a file; the program is exited with an error
   message if the file cannot be read or is invalid
*/
whitelist_t* whitelist_new(const char *path, Logger *logger);

/*
   Create a whitelist from <n> encoded barcodes of length <length> (see
   whitelist_encode())
*/
whitelist_t* whitelist_new_from_codes(const uint64_t *barcodes, size_t n,
    size_t length);

/*
   Encode the barcode <seq> (length: the barcode length); returns false if
   it contains characters other than ACGT
*/
bool whitelist_encode(const whitelist_t *self, const char *seq,
    uint64_t *barcode);

/* true if the encoded barcode is in the whitelist */
bool whitelist_contains(const whitelist_t *self, uint64_t barcode);

/* Delete a whitelist (NULL is accepted) */
void whitelist_delete(whitelist_t *self);

/* Length of the barcodes */
size_t whitelist_barcode_length(const whitelist_t *self);

/* Number of barcodes */
size_t whitelist_size(const whitelist_t *self);

/* Encoded barcodes, in no particular order (allocated, to be freed by the
   caller); returns their number */
size_t whitelist_codes(const whitelist_t *self, uint64_t **codes);

/*
   Search the barcode at the start of <seq> (length <len>) or, if
   <from_end>, at its end, with an edit distance of at most <max_dist>
   (<= WHITELIST_MAX_DIST); the result is stored in <match>.
*/
void whitelist_match(const whitelist_t *self, const char *seq, size_t len,
    bool from_end, long max_dist, whitelist_match_t *match);

/*
   Write the sequence of an encoded barcode to <out> (barcode length
   characters, plus a terminating null character)
*/
void whitelist_decode(const whitelist_t *self, uint64_t barcode, char *out);

/* Unit tests */
void whitelist_test(const char *tmpdir);

#endif
