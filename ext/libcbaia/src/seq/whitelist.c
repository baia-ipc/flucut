/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include "khash.h"
#include "alloc.h"
#include "logger.h"
#include "strutils.h"
#include "test.h"
#include "whitelist.h"

KHASH_SET_INIT_INT64(barcodes)

#define LINE_BUFSIZE 1024

struct whitelist_t {
  size_t length;
  khash_t(barcodes) *set;
};

static const char bases[] = "ACGT";

/* 2-bit code of a base, or -1 for other characters */
static inline int base_code(char c) {
  switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default: return -1;
  }
}

/* encode <len> characters; returns false if a character is not ACGT */
static bool encode(const char *seq, size_t len, uint64_t *code) {
  uint64_t c = 0;
  for (size_t i = 0; i < len; i++) {
    int b = base_code(seq[i]);
    if (b < 0)
      return false;
    c = (c << 2) | (uint64_t)b;
  }
  *code = c;
  return true;
}

void whitelist_decode(const whitelist_t *self, uint64_t barcode, char *out) {
  assert(self != NULL);
  assert(out != NULL);
  for (size_t i = self->length; i > 0; i--) {
    out[i - 1] = bases[barcode & 3];
    barcode >>= 2;
  }
  out[self->length] = '\0';
}

/* the set of the barcodes (exits if the allocation fails) */
static khash_t(barcodes) *new_set(void) {
  khash_t(barcodes) *set = kh_init(barcodes);
  if (set == NULL)
    fatal_printf("Failed allocating the whitelist\n");
  return set;
}

/* add <code> to <set> (exits if the allocation fails) */
static void add_barcode(khash_t(barcodes) *set, uint64_t code) {
  int ret;
  kh_put(barcodes, set, code, &ret);
  if (ret == -1)
    fatal_printf("Failed allocating the whitelist\n");
}

__attribute__ ((noreturn))
static void file_error(const char *path, size_t lineno, const char *msg,
    Logger *logger) {
  lprintf(logger, LOG_ERROR, "Whitelist %s, line %zu: %s\n", path, lineno,
          msg);
  exit(EXIT_FAILURE);
}

whitelist_t* whitelist_new(const char *path, Logger *logger) {
  assert(path != NULL);
  gzFile file = gzopen(path, "r");
  LCHECK(file != NULL, logger, "Opening file '%s' failed, %s\n", path,
         strerror(errno));
  whitelist_t *self;
  XMALLOC(self, sizeof(whitelist_t));
  self->length = 0;
  self->set = new_set();
  char line[LINE_BUFSIZE];
  size_t lineno = 0;
  while (gzgets(file, line, LINE_BUFSIZE) != NULL) {
    lineno++;
    size_t len = strlen(line);
    if (len == LINE_BUFSIZE - 1 && line[len - 1] != '\n')
      file_error(path, lineno, "line too long", logger);
    /* first word, without a "-" suffix */
    size_t from = 0;
    while (line[from] == ' ' || line[from] == '\t')
      from++;
    size_t to = from;
    while (line[to] != '\0' && line[to] != '\n' && line[to] != '\r' &&
           line[to] != ' ' && line[to] != '\t' && line[to] != '-')
      to++;
    if (to == from)
      continue;
    size_t blen = to - from;
    if (blen > WHITELIST_MAX_BARCODE_LENGTH) {
      char msg[64];
      snprintf(msg, sizeof(msg), "barcode longer than %d bases",
               WHITELIST_MAX_BARCODE_LENGTH);
      file_error(path, lineno, msg, logger);
    }
    if (self->length == 0)
      self->length = blen;
    else if (blen != self->length)
      file_error(path, lineno, "barcodes of different lengths",
                 logger);
    uint64_t code;
    if (!encode(line + from, blen, &code))
      file_error(path, lineno, "barcode with characters other than ACGT",
                 logger);
    add_barcode(self->set, code);
  }
  int errnum;
  const char *msg = gzerror(file, &errnum);
  LCHECK(errnum == Z_OK || errnum == Z_STREAM_END, logger,
         "Reading whitelist %s failed, %s\n", path,
         errnum == Z_ERRNO ? strerror(errno) : msg);
  gzclose(file);
  LCHECK(kh_size(self->set) > 0, logger, "Whitelist %s: no barcodes "
         "found\n", path);
  lprintf(logger, LOG_INFO, "Whitelist %s: %zu barcodes of length %zu\n",
          path, (size_t)kh_size(self->set), self->length);
  return self;
}

whitelist_t* whitelist_new_from_codes(const uint64_t *barcodes, size_t n,
    size_t length) {
  assert(barcodes != NULL || n == 0);
  assert(length > 0 && length <= WHITELIST_MAX_BARCODE_LENGTH);
  whitelist_t *self;
  XMALLOC(self, sizeof(whitelist_t));
  self->length = length;
  self->set = new_set();
  for (size_t i = 0; i < n; i++)
    add_barcode(self->set, barcodes[i]);
  return self;
}

bool whitelist_encode(const whitelist_t *self, const char *seq,
    uint64_t *barcode) {
  assert(self != NULL);
  assert(seq != NULL);
  assert(barcode != NULL);
  return encode(seq, self->length, barcode);
}

bool whitelist_contains(const whitelist_t *self, uint64_t barcode) {
  assert(self != NULL);
  return kh_get(barcodes, self->set, barcode) != kh_end(self->set);
}

void whitelist_delete(whitelist_t *self) {
  if (self == NULL)
    return;
  kh_destroy(barcodes, self->set);
  free(self);
}

size_t whitelist_barcode_length(const whitelist_t *self) {
  assert(self != NULL);
  return self->length;
}

size_t whitelist_size(const whitelist_t *self) {
  assert(self != NULL);
  return kh_size(self->set);
}

size_t whitelist_codes(const whitelist_t *self, uint64_t **codes) {
  assert(self != NULL);
  assert(codes != NULL);
  size_t n = 0;
  XMALLOC(*codes, sizeof(uint64_t) * (kh_size(self->set) + 1));
  for (khiter_t k = kh_begin(self->set); k != kh_end(self->set); k++)
    if (kh_exist(self->set, k))
      (*codes)[n++] = kh_key(self->set, k);
  return n;
}

/* --- search --- */

/* state of the enumeration of the sequences within a distance */
typedef struct {
  const whitelist_t *wl;
  const char *seq;
  size_t len;
  bool from_end;
  long max_cost;
  char buf[WHITELIST_MAX_BARCODE_LENGTH];
  /* results */
  size_t n_found;     /* distinct barcodes found (counted up to 2) */
  uint64_t barcode;
  size_t consumed;
} search_t;

/* i-th character of the sequence, from the anchored end */
static inline char seq_char(const search_t *s, size_t i) {
  return s->from_end ? s->seq[s->len - 1 - i] : s->seq[i];
}

/* position in the barcode of its j-th character from the anchored end */
static inline size_t buf_pos(const search_t *s, size_t j) {
  return s->from_end ? s->wl->length - 1 - j : j;
}

/* look up the barcode in the buffer, aligned to <consumed> characters */
static void lookup(search_t *s, size_t consumed) {
  uint64_t code;
  if (!encode(s->buf, s->wl->length, &code))
    return;
  if (kh_get(barcodes, s->wl->set, code) == kh_end(s->wl->set))
    return;
  if (s->n_found == 0) {
    s->n_found = 1;
    s->barcode = code;
    s->consumed = consumed;
  } else if (code == s->barcode) {
    /* same barcode, other alignment: prefer the aligned length closest
       to the barcode length, then the shortest */
    size_t l = s->wl->length;
    size_t d_new = consumed > l ? consumed - l : l - consumed,
           d_old = s->consumed > l ? s->consumed - l : l - s->consumed;
    if (d_new < d_old || (d_new == d_old && consumed < s->consumed))
      s->consumed = consumed;
  } else
    s->n_found = 2;
}

/* enumerate the barcodes obtained by editing the sequence: <i> characters
   of the sequence and <j> of the barcode are aligned, with cost <cost> */
static void enumerate(search_t *s, size_t i, size_t j, long cost) {
  size_t l = s->wl->length;
  if (j == l) {
    lookup(s, i);
    return;
  }
  if (cost == s->max_cost) {
    /* no further edits: copy the rest of the sequence */
    if (i + (l - j) > s->len)
      return;
    for (size_t k = 0; j + k < l; k++)
      s->buf[buf_pos(s, j + k)] = seq_char(s, i + k);
    lookup(s, i + (l - j));
    return;
  }
  size_t p = buf_pos(s, j);
  if (i < s->len) {
    /* match or substitution */
    char c = seq_char(s, i);
    for (size_t b = 0; b < 4; b++) {
      s->buf[p] = bases[b];
      enumerate(s, i + 1, j + 1, cost + (bases[b] == c ? 0 : 1));
    }
    /* character of the sequence not in the barcode */
    enumerate(s, i + 1, j, cost + 1);
  }
  /* character of the barcode not in the sequence */
  for (size_t b = 0; b < 4; b++) {
    s->buf[p] = bases[b];
    enumerate(s, i, j + 1, cost + 1);
  }
}

void whitelist_match(const whitelist_t *self, const char *seq, size_t len,
    bool from_end, long max_dist, whitelist_match_t *match) {
  assert(self != NULL);
  assert(seq != NULL || len == 0);
  assert(match != NULL);
  assert(max_dist >= 0 && max_dist <= WHITELIST_MAX_DIST);
  search_t s = {.wl = self, .seq = seq, .len = len, .from_end = from_end};
  match->status = WHITELIST_NOT_FOUND;
  /* the distances are tried in increasing order, thus all barcodes found
     for a distance are at exactly that distance */
  for (long d = 0; d <= max_dist; d++) {
    s.max_cost = d;
    s.n_found = 0;
    enumerate(&s, 0, 0, 0);
    if (s.n_found == 0)
      continue;
    if (s.n_found > 1)
      match->status = WHITELIST_AMBIGUOUS;
    else {
      match->status = d == 0 ? WHITELIST_EXACT : WHITELIST_CORRECTED;
      match->barcode = s.barcode;
      match->dist = d;
      match->consumed = s.consumed;
    }
    return;
  }
}

/* --- tests --- */

static const char *status_names[] = {"exact", "corrected", "ambiguous",
                                     "not_found"};

/* result as string: status, or barcode:dist:consumed */
static const char *match_str(const whitelist_t *wl, const char *seq,
    bool from_end, long max_dist) {
  static char out[64];
  whitelist_match_t m;
  whitelist_match(wl, seq, strlen(seq), from_end, max_dist, &m);
  if (m.status == WHITELIST_AMBIGUOUS || m.status == WHITELIST_NOT_FOUND)
    return status_names[m.status];
  char bc[WHITELIST_MAX_BARCODE_LENGTH + 1];
  whitelist_decode(wl, m.barcode, bc);
  snprintf(out, sizeof(out), "%s:%ld:%zu", bc, m.dist, m.consumed);
  return out;
}

void whitelist_test(const char *tmpdir) {
  char *path = strcat_alloc(tmpdir, "/cbaia-tests-whitelist.txt");
  FILE *f = fopen(path, "w");
  test_assert(f != NULL);
  fputs("AAAACCCCGGGGTTTT-1\nACGTACGTACGTACGT\n\nTTTTGGGGCCCCAAAA\n"
        "TTTTGGGGCCCCAAAC extra\n", f);
  fclose(f);
  whitelist_t *wl = whitelist_new(path, NULL);
  test_eq_size(whitelist_size(wl), 4);
  test_eq_size(whitelist_barcode_length(wl), 16);
  /* exact, at the start, followed by other characters */
  test_eq_str(match_str(wl, "AAAACCCCGGGGTTTTACGT", false, 2),
              "AAAACCCCGGGGTTTT:0:16");
  /* substitution */
  test_eq_str(match_str(wl, "AAAACCCCGAGGTTTTACGT", false, 2),
              "AAAACCCCGGGGTTTT:1:16");
  test_eq_str(match_str(wl, "AAAACCCCGAGGTTTTACGT", false, 0), "not_found");
  /* deletion in the sequence: one character less aligned */
  test_eq_str(match_str(wl, "AAAACCCCGGGTTTTACGT", false, 2),
              "AAAACCCCGGGGTTTT:1:15");
  /* insertion in the sequence: one character more aligned */
  test_eq_str(match_str(wl, "AAAACCCCAGGGGTTTTCGT", false, 2),
              "AAAACCCCGGGGTTTT:1:17");
  /* two errors */
  test_eq_str(match_str(wl, "AAAACCCAAGGGTTTTACG", false, 2),
              "AAAACCCCGGGGTTTT:2:16");
  test_eq_str(match_str(wl, "AAAACCCAAGGGTTTTACG", false, 1), "not_found");
  /* three errors (WHITELIST_MAX_DIST) */
  test_eq_str(match_str(wl, "AAAAGCCCGTGGTTCTACGT", false, 3),
              "AAAACCCCGGGGTTTT:3:16");
  test_eq_str(match_str(wl, "AAAAGCCCGTGGTTCTACGT", false, 2), "not_found");
  test_eq_str(match_str(wl, "CGTAAAAGCCCGTGGTTCT", true, 3),
              "AAAACCCCGGGGTTTT:3:16");
  /* N in the sequence */
  test_eq_str(match_str(wl, "ACGTACGTNCGTACGTAAA", false, 1),
              "ACGTACGTACGTACGT:1:16");
  /* two barcodes at distance 1 */
  test_eq_str(match_str(wl, "TTTTGGGGCCCCAAAG", false, 2), "ambiguous");
  /* at the end of the sequence */
  test_eq_str(match_str(wl, "GGGACGTACGTACGTACGT", true, 2),
              "ACGTACGTACGTACGT:0:16");
  test_eq_str(match_str(wl, "GGGACGTACGTACGTCGT", true, 2),
              "ACGTACGTACGTACGT:1:15");
  /* sequence shorter than the barcode */
  test_eq_str(match_str(wl, "ACGTACGTACGTACG", false, 1),
              "ACGTACGTACGTACGT:1:15");
  test_eq_str(match_str(wl, "", false, 2), "not_found");
  uint64_t codes[2];
  test_assert(whitelist_encode(wl, "ACGTACGTACGTACGT", codes));
  test_assert(whitelist_contains(wl, codes[0]));
  test_assert(!whitelist_encode(wl, "ACGTACGTACGTACGN", codes + 1));
  test_assert(whitelist_encode(wl, "ACGTACGTACGTACGG", codes + 1));
  test_assert(!whitelist_contains(wl, codes[1]));
  whitelist_t *wl2 = whitelist_new_from_codes(codes, 2, 16);
  test_eq_size(whitelist_size(wl2), 2);
  uint64_t *got;
  test_eq_size(whitelist_codes(wl2, &got), 2);
  test_assert(got[0] != got[1] && whitelist_contains(wl2, got[0]) &&
              whitelist_contains(wl2, got[1]));
  free(got);
  test_eq_str(match_str(wl2, "ACGTACGTACGTACGC", false, 1), "ambiguous");
  test_eq_str(match_str(wl2, "ACGTACGTACCTACGG", false, 1),
              "ACGTACGTACGTACGG:1:16");
  whitelist_delete(wl2);
  whitelist_delete(wl);
  unlink(path);
  free(path);
}
