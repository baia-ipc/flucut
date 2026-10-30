/* (c) 2026 Giorgio Gonnella */

/*
   MD5 message digest (RFC 1321), e.g. for the checksums of files in
   metadata (as the md5 of the onlists of seqspec files) or for checking
   downloads. Not to be used for security purposes.
*/

#ifndef MD5_H
#define MD5_H

#include <stddef.h>
#include <stdint.h>
#include "logger.h"

/* length of the hexadecimal digest, without the terminating null */
#define MD5_HEX_LEN 32

typedef struct {
  uint32_t state[4];
  uint64_t n_bytes;
  uint8_t buffer[64];
} md5_t;

/* Initialize the computation of a digest */
void md5_init(md5_t *self);

/* Add <len> bytes of <data> */
void md5_update(md5_t *self, const void *data, size_t len);

/* Finish the computation; write the digest as MD5_HEX_LEN lowercase
   hexadecimal digits, and a terminating null, to <hex> */
void md5_final(md5_t *self, char *hex);

/* Digest of the content of a file (to <hex>, as md5_final()); the program
   is exited with an error message if the file cannot be read */
void md5_file(const char *path, char *hex, Logger *logger);

/* Unit tests (<tmpdir>: directory for temporary files) */
void md5_test(const char *tmpdir);

#endif
