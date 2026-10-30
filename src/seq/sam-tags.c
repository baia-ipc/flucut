/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <ctype.h>
#include <string.h>
#include "sam-tags.h"
#include "test.h"

bool sam_tag_is_valid(const char *field, size_t len) {
  assert(field != NULL || len == 0);
  return len >= 5 && isalpha((unsigned char)field[0]) &&
         isalnum((unsigned char)field[1]) && field[2] == ':' &&
         strchr("AifZHB", field[3]) != NULL && field[4] == ':';
}

void sam_tag_start(strbuf_t *out, const char *tag, char type) {
  assert(out != NULL);
  assert(tag != NULL && strlen(tag) >= 2);
  char start[] = {'\t', tag[0], tag[1], ':', type, ':'};
  strbuf_append(out, start, sizeof(start));
}

void sam_tag_append_string(strbuf_t *out, const char *tag,
    const char *value, size_t len) {
  sam_tag_start(out, tag, 'Z');
  strbuf_append(out, value, len);
}

void sam_tag_append_int(strbuf_t *out, const char *tag, long value) {
  sam_tag_start(out, tag, 'i');
  strbuf_append_long(out, value);
}

void sam_tags_test(void) {
  test_assert(sam_tag_is_valid("RG:Z:run1", 9));
  test_assert(sam_tag_is_valid("ch:i:12", 7));
  test_assert(sam_tag_is_valid("M2:B:", 5));
  test_assert(!sam_tag_is_valid("ch:i:", 4));
  test_assert(!sam_tag_is_valid("ch=12", 5));
  test_assert(!sam_tag_is_valid("1h:i:12", 7));
  test_assert(!sam_tag_is_valid("ch:x:12", 7));
  test_assert(!sam_tag_is_valid("runid=abc", 9));
  strbuf_t out;
  strbuf_init(&out);
  sam_tag_append_string(&out, "CB", "ACGTX", 4);
  sam_tag_append_int(&out, "mk", 2);
  sam_tag_append_int(&out, "xn", -15);
  sam_tag_start(&out, "QT", 'Z');
  strbuf_append_str(&out, "II I");
  test_eq_mem(out.s, out.len, "\tCB:Z:ACGT\tmk:i:2\txn:i:-15\tQT:Z:II I");
  strbuf_free(&out);
}
