/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "alloc.h"
#include "bam-format.h"
#include "bam-reader.h"
#include "bgzf-compressor.h"
#include "sam-tags.h"
#include "strutils.h"
#include "test.h"
#include "ubam-writer.h"

void ubam_writer_header(strbuf_t *out, const char *text) {
  assert(out != NULL);
  assert(text != NULL);
  strbuf_append(out, BAM_MAGIC, BAM_MAGIC_LEN);
  bam_put32(out, (uint32_t)strlen(text));
  strbuf_append_str(out, text);
  bam_put32(out, 0);
}

/* 4-bit code of a base */
static uint8_t base_code(char c) {
  const char *p = strchr(BAM_BASES, toupper((unsigned char)c));
  return p == NULL || c == '\0' ? 15 : (uint8_t)(p - BAM_BASES);
}

/* parse an integer of a tag value (the whole text) */
static bool parse_int(const char *s, size_t len, long long *value) {
  char buf[32];
  if (len == 0 || len >= sizeof(buf))
    return false;
  memcpy(buf, s, len);
  buf[len] = '\0';
  char *end;
  errno = 0;
  *value = strtoll(buf, &end, 10);
  return errno == 0 && *end == '\0';
}

static bool parse_float(const char *s, size_t len, float *value) {
  char buf[64];
  if (len == 0 || len >= sizeof(buf))
    return false;
  memcpy(buf, s, len);
  buf[len] = '\0';
  char *end;
  errno = 0;
  *value = strtof(buf, &end);
  return errno == 0 && *end == '\0';
}

/* the smallest integer type for the value (0: out of range) */
static char int_type(long long v) {
  if (v < 0)
    return v >= INT8_MIN ? 'c' : v >= INT16_MIN ? 's' :
           v >= INT32_MIN ? 'i' : 0;
  return v <= UINT8_MAX ? 'C' : v <= UINT16_MAX ? 'S' :
         v <= UINT32_MAX ? 'I' : 0;
}

/* append a value of an integer or float type */
static void put_value(strbuf_t *out, char type, long long v, float f) {
  switch (type) {
    case 'c': case 'C': strbuf_append_char(out, (char)(v & 0xff)); break;
    case 's': case 'S': bam_put16(out, (uint16_t)(v & 0xffff)); break;
    case 'i': case 'I': bam_put32(out, (uint32_t)(v & 0xffffffff)); break;
    case 'f': {
      uint32_t u;
      memcpy(&u, &f, sizeof(u));
      bam_put32(out, u);
      break;
    }
    default:
      assert(false);
      break;
  }
}

/* append an array tag value (e.g. "C,1,2,3"); false if invalid (then
   <out> is unchanged) */
static bool put_array(strbuf_t *out, const char *v, size_t len) {
  if (len < 1 || strchr("cCsSiIf", v[0]) == NULL || v[0] == '\0')
    return false;
  char type = v[0];
  size_t start = out->len;
  strbuf_append_char(out, type);
  size_t count_pos = out->len;
  bam_put32(out, 0);
  uint32_t count = 0;
  const char *p = v + 1, *end = v + len;
  while (p < end) {
    if (*p != ',')
      goto invalid;
    p++;
    const char *e = memchr(p, ',', (size_t)(end - p));
    if (e == NULL)
      e = end;
    long long iv = 0;
    float fv = 0;
    if (type == 'f' ? !parse_float(p, (size_t)(e - p), &fv) :
        !parse_int(p, (size_t)(e - p), &iv))
      goto invalid;
    if (type != 'f') {
      long long lo = type == 'c' ? INT8_MIN : type == 's' ? INT16_MIN :
                     type == 'i' ? INT32_MIN : 0;
      long long hi = type == 'c' ? INT8_MAX : type == 'C' ? UINT8_MAX :
                     type == 's' ? INT16_MAX : type == 'S' ? UINT16_MAX :
                     type == 'i' ? INT32_MAX : UINT32_MAX;
      if (iv < lo || iv > hi)
        goto invalid;
    }
    put_value(out, type, iv, fv);
    count++;
    p = e;
  }
  bam_set32(out->s + count_pos, count);
  return true;
invalid:
  strbuf_truncate(out, start);
  return false;
}

/* size of the value of an integer or float type */
static size_t value_size(char type) {
  return type == 'c' || type == 'C' ? 1 : type == 's' || type == 'S' ? 2 : 4;
}

/* size of a tag written by put_tag() at <p> (name, type and value) */
static size_t tag_size(const char *p) {
  switch (p[2]) {
    case 'A': return 4;
    case 'Z': case 'H': return 3 + strlen(p + 3) + 1;
    case 'B': {
      const uint8_t *n = (const uint8_t *)p + 4;
      uint32_t count = (uint32_t)n[0] | ((uint32_t)n[1] << 8) |
                       ((uint32_t)n[2] << 16) | ((uint32_t)n[3] << 24);
      return 8 + value_size(p[3]) * count;
    }
    default: return 3 + value_size(p[2]);
  }
}

/* is there a tag with the name of <t> among the tags written by put_tag()
   in tags[0..len) */
static bool has_tag(const char *tags, size_t len, const char *t) {
  for (size_t i = 0; i < len; i += tag_size(tags + i))
    if (tags[i] == t[0] && tags[i + 1] == t[1])
      return true;
  return false;
}

/* append a tag given in SAM text format (e.g. "RG:Z:x"), if valid and not
   among the tags already written from out->s + <tags_from> */
static void put_tag(strbuf_t *out, const char *t, size_t len,
    size_t tags_from) {
  if (!sam_tag_is_valid(t, len) ||
      has_tag(out->s + tags_from, out->len - tags_from, t))
    return;
  const char *v = t + 5;
  size_t vlen = len - 5;
  size_t start = out->len;
  strbuf_append(out, t, 2);
  switch (t[3]) {
    case 'A':
      if (vlen != 1)
        goto invalid;
      strbuf_append_char(out, 'A');
      strbuf_append_char(out, v[0]);
      return;
    case 'i': {
      long long iv;
      char type;
      if (!parse_int(v, vlen, &iv) || (type = int_type(iv)) == 0)
        goto invalid;
      strbuf_append_char(out, type);
      put_value(out, type, iv, 0);
      return;
    }
    case 'f': {
      float fv;
      if (!parse_float(v, vlen, &fv))
        goto invalid;
      strbuf_append_char(out, 'f');
      put_value(out, 'f', 0, fv);
      return;
    }
    case 'Z': case 'H':
      if (memchr(v, '\0', vlen) != NULL)
        goto invalid;
      strbuf_append_char(out, t[3]);
      strbuf_append(out, v, vlen);
      strbuf_append_char(out, '\0');
      return;
    case 'B':
      strbuf_append_char(out, 'B');
      if (put_array(out, v, vlen))
        return;
      goto invalid;
    default:
      goto invalid;
  }
invalid:
  strbuf_truncate(out, start);
}

void ubam_writer_record(strbuf_t *out, const char *title, size_t title_len,
    const char *seq, const char *qual, size_t len) {
  assert(out != NULL);
  assert(title != NULL);
  /* the read name is the first word of the title (as in FASTQ files,
     whose titles are kept as they are, leading whitespace included) */
  const char *name = title, *end = title + title_len;
  while (name < end && is_space(*name))
    name++;
  size_t name_len = 0;
  while (name + name_len < end && !is_space(name[name_len]))
    name_len++;
  const char *tags = name + name_len;
  if (name_len == 0) {  /* (then there are no tags) */
    name = UBAM_NO_READ_NAME;
    name_len = strlen(UBAM_NO_READ_NAME);
  }
  size_t l_read_name = (name_len < UBAM_MAX_READ_NAME ? name_len :
                        UBAM_MAX_READ_NAME) + 1;
  size_t start = out->len;
  bam_put32(out, 0);   /* block_size, set at the end */
  bam_put32(out, 0xffffffff);
  bam_put32(out, 0xffffffff);
  strbuf_append_char(out, (char)l_read_name);
  strbuf_append_char(out, (char)255);
  bam_put16(out, BAM_UNMAPPED_BIN);
  bam_put16(out, 0);
  bam_put16(out, BAM_FUNMAP);
  bam_put32(out, (uint32_t)len);
  bam_put32(out, 0xffffffff);
  bam_put32(out, 0xffffffff);
  bam_put32(out, 0);
  strbuf_append(out, name, l_read_name - 1);
  strbuf_append_char(out, '\0');
  for (size_t i = 0; i < len; i += 2)
    strbuf_append_char(out, (char)((base_code(seq[i]) << 4) |
                       (i + 1 < len ? base_code(seq[i + 1]) : 0)));
  for (size_t i = 0; i < len; i++)
    strbuf_append_char(out, qual == NULL ? (char)BAM_QUAL_MISSING :
                       (char)bam_qual_value(qual[i]));
  const char *p = tags;
  size_t tags_from = out->len;
  while (p < end) {
    p++;
    const char *e = memchr(p, '\t', (size_t)(end - p));
    if (e == NULL)
      e = end;
    put_tag(out, p, (size_t)(e - p), tags_from);
    p = e;
  }
  uint32_t block_size = (uint32_t)(out->len - start - 4);
  bam_set32(out->s + start, block_size);
}

/* --- unit tests --- */

/* write a single record (no tags) to <b> (cleared before) */
static void write_record(strbuf_t *b, const char *title, const char *seq,
    const char *qual) {
  strbuf_clear(b);
  ubam_writer_record(b, title, strlen(title), seq, qual, strlen(seq));
}

/* the read name of the record in <b> (0-terminated) */
static const char *record_name(const strbuf_t *b) {
  return b->s + 4 + BAM_CORE_LEN;
}

/* the tags of the record in <b> (binary), and their length */
static const char *record_tags(const strbuf_t *b, size_t len,
    size_t *tags_len) {
  size_t from = 4 + BAM_CORE_LEN + (uint8_t)b->s[12] + (len + 1) / 2 + len;
  *tags_len = b->len - from;
  return b->s + from;
}

/* read name: the first word of the title */
static void test_read_name(void) {
  strbuf_t b;
  strbuf_init(&b);
  write_record(&b, " \tr1 x\tRG:Z:a", "A", "I");
  test_eq_str(record_name(&b), "r1");
  size_t tags_len;
  const char *tags = record_tags(&b, 1, &tags_len);
  test_assert(tags_len == 5 && memcmp(tags, "RGZa", 5) == 0);
  write_record(&b, "", "A", "I");
  test_eq_str(record_name(&b), UBAM_NO_READ_NAME);
  write_record(&b, "  ", "A", "I");
  test_eq_str(record_name(&b), UBAM_NO_READ_NAME);
  record_tags(&b, 1, &tags_len);
  test_eq_size(tags_len, 0);
  /* a name longer than UBAM_MAX_READ_NAME is truncated (the tags are
     kept) */
  char title[UBAM_MAX_READ_NAME + 20];
  memset(title, 'n', UBAM_MAX_READ_NAME + 10);
  strcpy(title + UBAM_MAX_READ_NAME + 10, "\tRG:Z:a");
  write_record(&b, title, "A", "I");
  test_eq_size(strlen(record_name(&b)), UBAM_MAX_READ_NAME);
  test_assert(strspn(record_name(&b), "n") == UBAM_MAX_READ_NAME);
  tags = record_tags(&b, 1, &tags_len);
  test_assert(tags_len == 5 && memcmp(tags, "RGZa", 5) == 0);
  /* exactly UBAM_MAX_READ_NAME characters: not truncated */
  strcpy(title + UBAM_MAX_READ_NAME, "\tRG:Z:a");
  write_record(&b, title, "A", "I");
  test_eq_size(strlen(record_name(&b)), UBAM_MAX_READ_NAME);
  strbuf_free(&b);
}

/* the qualities of the record in <b> (length <len>) */
static const uint8_t *record_qual(const strbuf_t *b, size_t len) {
  return (const uint8_t *)b->s + 4 + BAM_CORE_LEN + (uint8_t)b->s[12] +
         (len + 1) / 2;
}

/* qualities out of the range of SAM are clamped */
static void test_qualities(void) {
  strbuf_t b;
  strbuf_init(&b);
  write_record(&b, "r1", "ACGTA", " !I~\x80");
  const uint8_t *q = record_qual(&b, 5);
  test_assert(q[0] == 0 && q[1] == 0 && q[2] == 40 && q[3] == 93 &&
              q[4] == 93);
  /* missing qualities */
  strbuf_clear(&b);
  ubam_writer_record(&b, "r2", 2, "ACG", NULL, 3);
  q = record_qual(&b, 3);
  test_assert(q[0] == BAM_QUAL_MISSING && q[1] == BAM_QUAL_MISSING &&
              q[2] == BAM_QUAL_MISSING);
  strbuf_free(&b);
}

/* invalid and repeated tags are not written */
static void test_tags(void) {
  strbuf_t b;
  strbuf_init(&b);
  write_record(&b, "r1\tRG:Z:a\tXX:i:x\tRG:Z:b\tXX:i:2\tB1:B:C,1,2\t"
               "ch:i:5\tB1:B:c,3\tch:i:6\tXX:Z:y", "A", "I");
  size_t tags_len;
  const char *tags = record_tags(&b, 1, &tags_len);
  const char exp[] = "RGZa\0" "XXC\x02" "B1BC\x02\0\0\0\x01\x02" "chC\x05";
  test_eq_size(tags_len, sizeof(exp) - 1);
  test_assert(memcmp(tags, exp, sizeof(exp) - 1) == 0);
  /* a string value with a NUL character is not valid */
  strbuf_clear(&b);
  ubam_writer_record(&b, "r1\tst:Z:a\0b", 11, "A", "I", 1);
  record_tags(&b, 1, &tags_len);
  test_eq_size(tags_len, 0);
  /* a name longer than UBAM_MAX_READ_NAME is truncated (the tags are
     kept) */
  char title[UBAM_MAX_READ_NAME + 20];
  memset(title, 'n', UBAM_MAX_READ_NAME + 10);
  strcpy(title + UBAM_MAX_READ_NAME + 10, "\tRG:Z:a");
  write_record(&b, title, "A", "I");
  test_eq_size(strlen(record_name(&b)), UBAM_MAX_READ_NAME);
  test_assert(strspn(record_name(&b), "n") == UBAM_MAX_READ_NAME);
  tags = record_tags(&b, 1, &tags_len);
  test_assert(tags_len == 5 && memcmp(tags, "RGZa", 5) == 0);
  /* exactly UBAM_MAX_READ_NAME characters: not truncated */
  strcpy(title + UBAM_MAX_READ_NAME, "\tRG:Z:a");
  write_record(&b, title, "A", "I");
  test_eq_size(strlen(record_name(&b)), UBAM_MAX_READ_NAME);
  strbuf_free(&b);
}

void ubam_writer_test(const char *tmpdir) {
  test_read_name();
  test_qualities();
  test_tags();
  strbuf_t data, out;
  strbuf_init(&data);
  strbuf_init(&out);
  ubam_writer_header(&data, "@HD\tVN:1.6\tSO:unknown\n");
  const char *title = "r1\tqs:i:12\tch:i:-300\tdu:f:1.5\tst:Z:a b\t"
                      "XA:A:x\tbad\tXB:B:c,-1,2\tXI:i:x\tXC:i:4294967295\t"
                      "CB:Z:ACGT";
  ubam_writer_record(&data, title, strlen(title), "acgtn", "ABCDE", 5);
  ubam_writer_record(&data, "r2", 2, "ACG", "!!I", 3);
  ubam_writer_record(&data, "r3", 2, "ACG", NULL, 3);
  bgzf_compressor_t *c = bgzf_compressor_new(6);
  bgzf_compressor_run(c, data.s, data.len, &out);
  bgzf_compressor_delete(c);
  strbuf_append(&out, (const char *)bgzf_eof, BGZF_EOF_LEN);
  char *path = strcat_alloc(tmpdir, "/cbaia-tests-out.bam");
  FILE *f = fopen(path, "wb");
  test_assert(f != NULL);
  fwrite(out.s, 1, out.len, f);
  fclose(f);
  test_assert(bam_reader_has_eof_marker(path));
  bam_reader_t *r = bam_reader_new(path, NULL);
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  test_eq(bam_reader_read_batch(r, &batch, 10), (size_t)3, size_t, "%zu");
  /* invalid fields are not written, the array tag is not read */
  const char *exp = "r1\tqs:i:12\tch:i:-300\tdu:f:1.5\tst:Z:a b\tXA:A:x\t"
                    "XC:i:4294967295\tCB:Z:ACGT";
  const fastq_record_t *rec = batch.records;
  test_eq(rec->title_len, strlen(exp), size_t, "%zu");
  test_assert(memcmp(batch.data.s + rec->title_from, exp, strlen(exp)) == 0);
  test_assert(memcmp(batch.data.s + rec->seq_from, "ACGTN", 5) == 0);
  test_assert(memcmp(batch.data.s + rec->qual_from, "ABCDE", 5) == 0);
  test_assert(!rec->qual_missing);
  /* qualities 0 are not missing qualities */
  rec = batch.records + 1;
  test_eq(rec->title_len, (size_t)2, size_t, "%zu");
  test_assert(memcmp(batch.data.s + rec->qual_from, "!!I", 3) == 0);
  test_assert(!rec->qual_missing);
  /* missing qualities: read as '!', marked, and written again as missing */
  rec = batch.records + 2;
  test_assert(memcmp(batch.data.s + rec->qual_from, "!!!", 3) == 0);
  test_assert(rec->qual_missing);
  strbuf_t again;
  strbuf_init(&again);
  ubam_writer_record(&again, batch.data.s + rec->title_from, rec->title_len,
                     batch.data.s + rec->seq_from, rec->qual_missing ? NULL :
                     batch.data.s + rec->qual_from, rec->seq_len);
  const uint8_t *q = record_qual(&again, 3);
  test_assert(q[0] == BAM_QUAL_MISSING && q[1] == BAM_QUAL_MISSING &&
              q[2] == BAM_QUAL_MISSING);
  strbuf_free(&again);
  fastq_batch_free(&batch);
  bam_reader_delete(r);
  unlink(path);
  free(path);
  strbuf_free(&data);
  strbuf_free(&out);
}
