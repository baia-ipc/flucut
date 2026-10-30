/* (c) 2026 Giorgio Gonnella */

/*
   Encoding of unaligned BAM data (header and records), to be compressed
   with BGZF (see bgzf-compressor.h).

   A record is given as a FASTQ-like entry, whose title contains the read
   name and SAM tags in text format, separated by tabs (as in the FASTQ
   output of Dorado, or of bam-reader.h): the record is unmapped (flag 4),
   the tags are encoded in binary form (integers with the smallest type,
   as in htslib); fields of the title which are not valid SAM tags are not
   written, nor the tags with the name of a tag already written (a record
   cannot have two tags with the same name: the first valid one is kept).
*/

#ifndef UBAM_WRITER_H
#define UBAM_WRITER_H

#include <stddef.h>
#include "strbuf.h"

/* maximum length of a read name in BAM (longer names are truncated) */
#define UBAM_MAX_READ_NAME 254

/* read name written for an empty title (the unavailable QNAME of SAM) */
#define UBAM_NO_READ_NAME "*"

/* Append the header (SAM header text <text>, no reference) to <out> */
void ubam_writer_header(strbuf_t *out, const char *text);

/*
   Append a record to <out>: <title> (length <title_len>): read name (the
   first word, leading whitespace is skipped; UBAM_NO_READ_NAME if none)
   and tags; <seq> and <qual> (Phred+33; characters below '!' or above
   '~' are clamped to these; NULL: missing qualities, written as 0xFF):
   <len> bases. The BAM reader gives missing qualities as '!' and marks
   the record (qual_missing, see bam-reader.h): for such a record, pass
   NULL to keep them missing.
*/
void ubam_writer_record(strbuf_t *out, const char *title, size_t title_len,
    const char *seq, const char *qual, size_t len);

/* Unit tests (<tmpdir>: directory for temporary files) */
void ubam_writer_test(const char *tmpdir);

#endif
