/* (c) 2026 Giorgio Gonnella */

/*
   Writing of FASTQ records to a buffer (e.g. for a batch of records,
   which is then compressed with gzip-compressor.h and written to a file).
*/

#ifndef FASTQ_WRITER_H
#define FASTQ_WRITER_H

#include <stddef.h>
#include "strbuf.h"

/*
   Append a record to <out>: the title <title> (length <title_len>,
   without "@"), the sequence <seq> and the qualities <qual> (<len>
   characters each), in four lines
*/
void fastq_writer_record(strbuf_t *out, const char *title, size_t title_len,
    const char *seq, const char *qual, size_t len);

/* Unit tests */
void fastq_writer_test(void);

#endif
