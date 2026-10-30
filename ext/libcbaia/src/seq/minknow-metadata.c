/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <ctype.h>
#include <string.h>
#include "minknow-metadata.h"
#include "sam-tags.h"
#include "test.h"

/* the read metadata of MinKNOW: key, SAM tag, integer value */
static const struct {
  const char *key, *tag;
  bool integer;
} minknow_fields[] = {
  {"runid", "ri", false}, {"read", "rn", true}, {"ch", "ch", true},
  {"start_time", "st", false}, {"flow_cell_id", "PU", false},
  {"protocol_group_id", "pg", false}, {"sample_id", "si", false},
  {"barcode", MINKNOW_BARCODE_SAM_TAG, false}, {"barcode_alias", "oa", false},
  {"parent_read_id", "pi", false},
  {"basecall_model_version_id", "bm", false}};
#define N_MINKNOW_FIELDS (sizeof(minknow_fields) / sizeof(minknow_fields[0]))

static bool all_digits(const char *s, size_t len) {
  if (len == 0)
    return false;
  for (size_t i = 0; i < len; i++)
    if (!isdigit((unsigned char)s[i]))
      return false;
  return true;
}

bool minknow_metadata_append_sam_tag(strbuf_t *out, const char *field,
    size_t len) {
  assert(out != NULL);
  assert(field != NULL || len == 0);
  const char *eq = memchr(field, '=', len);
  if (eq == NULL)
    return false;
  size_t key_len = (size_t)(eq - field), value_len = len - key_len - 1;
  for (size_t i = 0; i < N_MINKNOW_FIELDS; i++)
    if (strlen(minknow_fields[i].key) == key_len &&
        memcmp(field, minknow_fields[i].key, key_len) == 0) {
      bool integer = minknow_fields[i].integer &&
                     all_digits(eq + 1, value_len);
      sam_tag_start(out, minknow_fields[i].tag, integer ? 'i' : 'Z');
      strbuf_append(out, eq + 1, value_len);
      return true;
    }
  return false;
}

/* conversion of a field (0-terminated), see
   minknow_metadata_append_sam_tag */
static bool convert(strbuf_t *out, const char *field) {
  return minknow_metadata_append_sam_tag(out, field, strlen(field));
}

void minknow_metadata_test(void) {
  strbuf_t out;
  strbuf_init(&out);
  test_assert(convert(&out, "runid=abc"));
  test_assert(convert(&out, "read=12"));
  test_assert(convert(&out, "ch=7"));
  test_assert(convert(&out, "start_time=2024-05-01T10:00:00Z"));
  test_assert(convert(&out, "flow_cell_id=FAX12345"));
  test_assert(convert(&out, "protocol_group_id=exp1"));
  test_assert(convert(&out, "sample_id=pool1"));
  test_assert(convert(&out, "barcode=barcode01"));
  test_assert(convert(&out, "barcode_alias=lib1"));
  test_assert(convert(&out, "parent_read_id=r1"));
  test_assert(convert(&out, "basecall_model_version_id=dna_r10@v5"));
  test_eq_mem(out.s, out.len, "\tri:Z:abc\trn:i:12\tch:i:7"
      "\tst:Z:2024-05-01T10:00:00Z\tPU:Z:FAX12345\tpg:Z:exp1\tsi:Z:pool1"
      "\tob:Z:barcode01\toa:Z:lib1\tpi:Z:r1\tbm:Z:dna_r10@v5");
  /* not converted: unknown key, no '=', key prefix of a known one */
  strbuf_clear(&out);
  test_assert(!convert(&out, "newkey=1"));
  test_assert(!convert(&out, "word"));
  test_assert(!convert(&out, "run=abc"));
  test_assert(!convert(&out, "runidx=abc"));
  test_eq_size(out.len, 0);
  /* integer fields with another value, empty value */
  test_assert(convert(&out, "ch=x1"));
  test_assert(convert(&out, "read="));
  test_assert(convert(&out, "sample_id="));
  test_eq_mem(out.s, out.len, "\tch:Z:x1\trn:Z:\tsi:Z:");
  strbuf_free(&out);
}
