/* (c) 2026 Giorgio Gonnella */

/*
   BGZF compression (the format of BAM files): the data are compressed in
   blocks of at most 64 KiB, each a gzip member with an extra field (BC)
   containing the size of the block. As for the gzip members of
   gzip-compressor.h, blocks can be compressed in parallel and written in
   order; a BGZF file ends with an empty block (the end-of-file marker).
*/

#ifndef BGZF_COMPRESSOR_H
#define BGZF_COMPRESSOR_H

#include <stddef.h>
#include <stdint.h>
#include "strbuf.h"

/* the end-of-file marker: an empty BGZF block */
#define BGZF_EOF_LEN 28
extern const uint8_t bgzf_eof[BGZF_EOF_LEN];

typedef struct bgzf_compressor_t bgzf_compressor_t;

/* Create a new compressor (use one for each thread) */
bgzf_compressor_t* bgzf_compressor_new(int level);

/* Compress <len> bytes of <data> as BGZF blocks, appended to <out> */
void bgzf_compressor_run(bgzf_compressor_t *self, const char *data,
    size_t len, strbuf_t *out);

/* Delete a compressor (NULL is accepted) */
void bgzf_compressor_delete(bgzf_compressor_t *self);

/* Unit tests (<tmpdir>: directory for temporary files) */
void bgzf_compressor_test(const char *tmpdir);

#endif
