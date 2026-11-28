/* (c) 2026 Giorgio Gonnella */

/*
   BAM data and files for the unit tests of the modules reading or writing
   BAM files (internal: not part of the API, it can change).
*/

#ifndef TEST_BAM_H
#define TEST_BAM_H

#include <stdbool.h>
#include <stddef.h>
#include "strbuf.h"

/*
   Append to <out> BGZF-compressed unaligned BAM data: the header with the
   SAM header text <text>, if not NULL, then a record for each of the
   <n_titles> <titles> (sequence ACGT, qualities IIII)
*/
void test_bam_data(strbuf_t *out, const char *text,
    const char *const *titles, size_t n_titles);

/*
   Write the BAM file <tmpdir>/<name>, with the header <text> and the
   records of test_bam_data(), followed by the BGZF end-of-file marker if
   <eof>; returns the path (to be freed by the caller)
*/
char *test_bam_file(const char *tmpdir, const char *name, const char *text,
    const char *const *titles, size_t n_titles, bool eof);

#endif
