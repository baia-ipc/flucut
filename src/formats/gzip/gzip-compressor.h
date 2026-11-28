/* (c) 2026 Giorgio Gonnella */

/*
   Compression of data blocks as gzip members; the concatenation of gzip
   members is a valid gzip file, thus blocks can be compressed in parallel
   and written in order.
*/

#ifndef GZIP_COMPRESSOR_H
#define GZIP_COMPRESSOR_H

#include <stddef.h>
#include "strbuf.h"

typedef struct gzip_compressor_t gzip_compressor_t;

/* Create a new compressor (use one for each thread) */
gzip_compressor_t* gzip_compressor_new(int level);

/*
   Compress <len> bytes of <data> as a gzip member, appended to <out>
*/
void gzip_compressor_run(gzip_compressor_t *self, const char *data,
    size_t len, strbuf_t *out);

/* Delete a compressor (NULL is accepted) */
void gzip_compressor_delete(gzip_compressor_t *self);

/* Unit tests (<tmpdir>: directory for temporary files) */
void gzip_compressor_test(const char *tmpdir);

#endif
