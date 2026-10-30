/* (c) 2026 Giorgio Gonnella */

/*
   Implementation of RFC 1321: the message is processed in blocks of 64
   bytes, each in 4 rounds of 16 operations; it is padded with a 1 bit,
   zeros and its length in bits (little endian).
*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "alloc.h"
#include "file.h"
#include "md5.h"
#include "strutils.h"
#include "test.h"

#define MD5_BUFSIZE (1 << 16)

/* the shift amounts of each operation */
static const uint32_t s_shift[64] = {
  7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22,
  5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20,
  4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23,
  6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21,
};

/* the constants: floor(abs(sin(i + 1)) * 2^32) */
static const uint32_t k_const[64] = {
  0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee, 0xf57c0faf, 0x4787c62a,
  0xa8304613, 0xfd469501, 0x698098d8, 0x8b44f7af, 0xffff5bb1, 0x895cd7be,
  0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821, 0xf61e2562, 0xc040b340,
  0x265e5a51, 0xe9b6c7aa, 0xd62f105d, 0x02441453, 0xd8a1e681, 0xe7d3fbc8,
  0x21e1cde6, 0xc33707d6, 0xf4d50d87, 0x455a14ed, 0xa9e3e905, 0xfcefa3f8,
  0x676f02d9, 0x8d2a4c8a, 0xfffa3942, 0x8771f681, 0x6d9d6122, 0xfde5380c,
  0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70, 0x289b7ec6, 0xeaa127fa,
  0xd4ef3085, 0x04881d05, 0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665,
  0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039, 0x655b59c3, 0x8f0ccc92,
  0xffeff47d, 0x85845dd1, 0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1,
  0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391,
};

static uint32_t rotl(uint32_t x, uint32_t c) {
  return (x << c) | (x >> (32 - c));
}

/* process a block of 64 bytes */
static void block(uint32_t state[4], const uint8_t *p) {
  uint32_t m[16];
  for (size_t i = 0; i < 16; i++)
    m[i] = (uint32_t)p[4 * i] | ((uint32_t)p[4 * i + 1] << 8) |
           ((uint32_t)p[4 * i + 2] << 16) | ((uint32_t)p[4 * i + 3] << 24);
  uint32_t a = state[0], b = state[1], c = state[2], d = state[3];
  for (uint32_t i = 0; i < 64; i++) {
    uint32_t f, g;
    if (i < 16) {
      f = (b & c) | (~b & d);
      g = i;
    } else if (i < 32) {
      f = (d & b) | (~d & c);
      g = (5 * i + 1) % 16;
    } else if (i < 48) {
      f = b ^ c ^ d;
      g = (3 * i + 5) % 16;
    } else {
      f = c ^ (b | ~d);
      g = (7 * i) % 16;
    }
    uint32_t tmp = d;
    d = c;
    c = b;
    b = b + rotl(a + f + k_const[i] + m[g], s_shift[i]);
    a = tmp;
  }
  state[0] += a;
  state[1] += b;
  state[2] += c;
  state[3] += d;
}

void md5_init(md5_t *self) {
  assert(self != NULL);
  self->state[0] = 0x67452301;
  self->state[1] = 0xefcdab89;
  self->state[2] = 0x98badcfe;
  self->state[3] = 0x10325476;
  self->n_bytes = 0;
}

void md5_update(md5_t *self, const void *data, size_t len) {
  assert(self != NULL);
  assert(data != NULL || len == 0);
  const uint8_t *p = data;
  size_t used = (size_t)(self->n_bytes % 64);
  self->n_bytes += len;
  if (used > 0) {
    size_t n = 64 - used < len ? 64 - used : len;
    memcpy(self->buffer + used, p, n);
    p += n;
    len -= n;
    if (used + n < 64)
      return;
    block(self->state, self->buffer);
  }
  for (; len >= 64; p += 64, len -= 64)
    block(self->state, p);
  memcpy(self->buffer, p, len);
}

void md5_final(md5_t *self, char *hex) {
  assert(self != NULL);
  assert(hex != NULL);
  uint64_t bits = self->n_bytes * 8;
  uint8_t pad[72] = {0x80};
  size_t used = (size_t)(self->n_bytes % 64);
  size_t pad_len = used < 56 ? 56 - used : 120 - used;
  for (size_t i = 0; i < 8; i++)
    pad[pad_len + i] = (uint8_t)(bits >> (8 * i));
  md5_update(self, pad, pad_len + 8);
  static const char digits[] = "0123456789abcdef";
  for (size_t i = 0; i < 16; i++) {
    uint8_t byte = (uint8_t)(self->state[i / 4] >> (8 * (i % 4)));
    hex[2 * i] = digits[byte >> 4];
    hex[2 * i + 1] = digits[byte & 15];
  }
  hex[MD5_HEX_LEN] = '\0';
}

void md5_file(const char *path, char *hex, Logger *logger) {
  assert(path != NULL);
  FILE *f = xfopen(path, NULL, "rb", "file (MD5 digest)", logger);
  md5_t md5;
  md5_init(&md5);
  uint8_t *buf;
  XMALLOC(buf, MD5_BUFSIZE);
  size_t n;
  while ((n = fread(buf, 1, MD5_BUFSIZE, f)) > 0)
    md5_update(&md5, buf, n);
  LCHECK(!ferror(f), logger, "Reading file '%s' failed, %s\n", path,
         strerror(errno));
  free(buf);
  xfclose(&f, logger);
  md5_final(&md5, hex);
}

/* --- unit tests --- */

static void check(const char *s, const char *expected) {
  md5_t md5;
  char hex[MD5_HEX_LEN + 1];
  md5_init(&md5);
  md5_update(&md5, s, strlen(s));
  md5_final(&md5, hex);
  test_eq_str(hex, expected);
  /* the same, added one byte at a time */
  md5_init(&md5);
  for (size_t i = 0; s[i] != '\0'; i++)
    md5_update(&md5, s + i, 1);
  md5_final(&md5, hex);
  test_eq_str(hex, expected);
}

void md5_test(const char *tmpdir) {
  /* the test suite of RFC 1321 */
  check("", "d41d8cd98f00b204e9800998ecf8427e");
  check("a", "0cc175b9c0f1b6a831c399e269772661");
  check("abc", "900150983cd24fb0d6963f7d28e17f72");
  check("message digest", "f96b697d7cb7938d525a2f31aaf161d0");
  check("abcdefghijklmnopqrstuvwxyz", "c3fcd3d76192e4007dfb496cca67e13b");
  check("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789",
        "d174ab98d277d9f5a5611c2c9f419d9f");
  check("1234567890123456789012345678901234567890"
        "1234567890123456789012345678901234567890",
        "57edf4a22be3c955ac49da2e2107b67a");
  /* a file larger than the buffer: 200,000 times "abc" */
  char *path = strcat_alloc(tmpdir, "/cbaia-tests-md5.txt");
  FILE *f = fopen(path, "w");
  test_assert(f != NULL);
  md5_t md5;
  md5_init(&md5);
  for (size_t i = 0; i < 200000; i++) {
    fputs("abc", f);
    md5_update(&md5, "abc", 3);
  }
  fclose(f);
  char expected[MD5_HEX_LEN + 1], hex[MD5_HEX_LEN + 1];
  md5_final(&md5, expected);
  md5_file(path, hex, NULL);
  test_eq_str(hex, expected);
  unlink(path);
  free(path);
}
