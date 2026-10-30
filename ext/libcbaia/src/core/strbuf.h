/* (c) 2026 Giorgio Gonnella */

/*
   Growable string buffer, for building output text without allocations
   for each piece (the memory is kept when the buffer is cleared).
*/

#ifndef STRBUF_H
#define STRBUF_H

#include <assert.h>
#include <stddef.h>
#include <string.h>
#include "alloc.h"

/* minimum free space after a reallocation; the size is also increased
   by 50%, for amortized constant time appends */
#define STRBUF_INC 4096

typedef struct {
  char *s;     /* not 0-terminated */
  size_t len;
  size_t cap;
} strbuf_t;

/* Initialize an empty buffer */
static inline void strbuf_init(strbuf_t *self) {
  self->s = NULL;
  self->len = 0;
  self->cap = 0;
}

/* Free the memory of the buffer */
static inline void strbuf_free(strbuf_t *self) {
  free(self->s);
  strbuf_init(self);
}

/* Remove the content (keeping the memory) */
static inline void strbuf_clear(strbuf_t *self) {
  self->len = 0;
}

/* Remove the content after the first <len> characters (at most the
   length of the content), keeping the memory */
static inline void strbuf_truncate(strbuf_t *self, size_t len) {
  assert(len <= self->len);
  self->len = len;
}

/* Append <len> characters */
static inline void strbuf_append(strbuf_t *self, const char *s,
    size_t len) {
  if (len == 0)  /* also if <s> is NULL */
    return;
  XCHECKSPACEINCPERC(self->s, &self->cap, self->len + len + STRBUF_INC, 50);
  memcpy(self->s + self->len, s, len);
  self->len += len;
}

/* Append a 0-terminated string */
static inline void strbuf_append_str(strbuf_t *self, const char *s) {
  strbuf_append(self, s, strlen(s));
}

/* Append a character */
static inline void strbuf_append_char(strbuf_t *self, char c) {
  XCHECKSPACEINCPERC(self->s, &self->cap, self->len + 1 + STRBUF_INC, 50);
  self->s[self->len++] = c;
}

/*
   Space for writing up to <n> characters directly at the end of the
   content (e.g. the output of a compressor): returns a pointer to it; the
   characters written are then added to the content by strbuf_commit()
*/
static inline char *strbuf_reserve(strbuf_t *self, size_t n) {
  XCHECKSPACEINCPERC(self->s, &self->cap, self->len + n, 50);
  return self->s + self->len;
}

/* Add to the content <n> characters written in the space given by
   strbuf_reserve() (at most those reserved) */
static inline void strbuf_commit(strbuf_t *self, size_t n) {
  self->len += n;
}

/* Append an integer in decimal notation */
static inline void strbuf_append_long(strbuf_t *self, long value) {
  char digits[24];
  size_t n = 0;
  unsigned long v = value < 0 ? -(unsigned long)value : (unsigned long)value;
  do {
    digits[n++] = (char)('0' + v % 10);
    v /= 10;
  } while (v > 0);
  if (value < 0)
    digits[n++] = '-';
  XCHECKSPACEINCPERC(self->s, &self->cap, self->len + n + STRBUF_INC, 50);
  while (n > 0)
    self->s[self->len++] = digits[--n];
}

#endif
