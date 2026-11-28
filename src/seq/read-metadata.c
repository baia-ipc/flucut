/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <string.h>
#include "read-metadata.h"
#include "test.h"

static const char *const keys[READ_METADATA_N] = {
  "flow_cell", "run", "basecall_model", "sample", "experiment",
  "instrument", "start_time"};
static const char *const labels[READ_METADATA_N] = {
  "flow cell", "run", "basecalling model", "sample name", "experiment",
  "instrument", "first read"};
/* the fields of the titles written by MinKNOW (NULL: none) */
static const char *const fields[READ_METADATA_N] = {
  "flow_cell_id=", "runid=", "basecall_model_version_id=", "sample_id=",
  "protocol_group_id=", NULL, "start_time="};

const char* read_metadata_key(size_t i) {
  assert(i < READ_METADATA_N);
  return keys[i];
}

const char* read_metadata_label(size_t i) {
  assert(i < READ_METADATA_N);
  return labels[i];
}

void read_metadata_clear(read_metadata_t *self) {
  assert(self != NULL);
  memset(self, 0, sizeof(read_metadata_t));
}

/* add the value <value>, of length <len>, of the value <i> */
static void add_value(read_metadata_t *self, size_t i, const char *value,
    size_t len) {
  if (len == 0)
    return;
  if (len >= READ_METADATA_VALUE_SIZE)
    len = READ_METADATA_VALUE_SIZE - 1;
  char *known = self->values[i];
  bool same = strlen(known) == len && memcmp(known, value, len) == 0;
  if (known[0] != '\0' && same)
    return;
  /* (the times, in the same format, are ordered as their text) */
  if (known[0] != '\0' && i == READ_METADATA_START &&
      strncmp(known, value, len) <= 0)
    return;
  if (known[0] != '\0' && i != READ_METADATA_START)
    self->several[i] = true;
  memcpy(known, value, len);
  known[len] = '\0';
}

/* the length of the field starting at <s>, which ends before <end>: up
   to a space or a tab */
static size_t field_length(const char *s, const char *end) {
  size_t len = 0;
  while (s + len < end && s[len] != ' ' && s[len] != '\t')
    len++;
  return len;
}

/* add the values of the first word of a title of an Illumina read
   (<instrument>:<run>:<flow cell>:<lane>:<tile>:<x>:<y>), if it is one */
static void add_illumina(read_metadata_t *self, const char *word,
    size_t len) {
  const char *parts[3];
  size_t lengths[3], n = 0, from = 0;
  for (size_t i = 0; i <= len; i++) {
    if (i < len && word[i] != ':')
      continue;
    if (n < 3) {
      parts[n] = word + from;
      lengths[n] = i - from;
    }
    n++;
    from = i + 1;
  }
  if (n != 7)
    return;
  add_value(self, READ_METADATA_INSTRUMENT, parts[0], lengths[0]);
  add_value(self, READ_METADATA_RUN, parts[1], lengths[1]);
  add_value(self, READ_METADATA_FLOW_CELL, parts[2], lengths[2]);
}

void read_metadata_add_title(read_metadata_t *self, const char *title,
    size_t len) {
  assert(self != NULL);
  assert(title != NULL);
  const char *end = title + len;
  if (len > 0 && title[0] == '@')
    title++;
  size_t word = field_length(title, end);
  bool any = false;
  for (const char *s = title + word + 1; s < end; s++) {
    /* (a field starts after a space or a tab; the first one can start
       after the separator which ends the first word) */
    if (s[-1] != ' ' && s[-1] != '\t')
      continue;
    for (size_t i = 0; i < READ_METADATA_N; i++) {
      size_t n = fields[i] == NULL ? 0 : strlen(fields[i]);
      if (n == 0 || (size_t)(end - s) < n || memcmp(s, fields[i], n) != 0)
        continue;
      add_value(self, i, s + n, field_length(s + n, end));
      any = true;
    }
  }
  if (!any)
    add_illumina(self, title, word);
}

void read_metadata_merge(read_metadata_t *self,
    const read_metadata_t *other) {
  assert(self != NULL);
  assert(other != NULL);
  for (size_t i = 0; i < READ_METADATA_N; i++) {
    add_value(self, i, other->values[i], strlen(other->values[i]));
    self->several[i] = self->several[i] || other->several[i];
  }
}

bool read_metadata_known(const read_metadata_t *self) {
  assert(self != NULL);
  for (size_t i = 0; i < READ_METADATA_N; i++)
    if (self->values[i][0] != '\0')
      return true;
  return false;
}

const char* read_metadata_text(const read_metadata_t *self, size_t i) {
  assert(self != NULL);
  assert(i < READ_METADATA_N);
  return self->several[i] ? "several" :
         self->values[i][0] != '\0' ? self->values[i] : NULL;
}

void read_metadata_save(const read_metadata_t *self, FILE *file) {
  assert(self != NULL);
  assert(file != NULL);
  for (size_t i = 0; i < READ_METADATA_N; i++)
    if (self->values[i][0] != '\0')
      fprintf(file, "meta %zu %d %s\n", i, self->several[i] ? 1 : 0,
              self->values[i]);
}

bool read_metadata_load_line(read_metadata_t *self, const char *line) {
  assert(self != NULL);
  assert(line != NULL);
  size_t i;
  int several, n;
  if (sscanf(line, "meta %zu %d %n", &i, &several, &n) != 2 ||
      i >= READ_METADATA_N)
    return false;
  const char *value = line + n;
  add_value(self, i, value, strcspn(value, "\n"));
  self->several[i] = self->several[i] || several != 0;
  return true;
}

void read_metadata_test(void) {
  read_metadata_t m, other;
  read_metadata_clear(&m);
  test_assert(!read_metadata_known(&m));
  const char *t1 = "@r1 runid=abc ch=143 "
    "start_time=2025-05-28T11:41:01.294251+07:00 flow_cell_id=FBB01520 "
    "sample_id=flu12 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_sup";
  read_metadata_add_title(&m, t1, strlen(t1));
  test_assert(read_metadata_known(&m));
  test_eq_str(read_metadata_text(&m, READ_METADATA_FLOW_CELL), "FBB01520");
  test_eq_str(read_metadata_text(&m, READ_METADATA_RUN), "abc");
  test_eq_str(read_metadata_text(&m, READ_METADATA_MODEL),
              "dna_r10.4.1_e8.2_400bps_sup");
  test_assert(read_metadata_text(&m, READ_METADATA_EXPERIMENT) == NULL);
  /* a read of another run, earlier: the run is then several, the
     earliest start is kept */
  const char *t2 = "r2 runid=def start_time=2025-05-28T10:00:00+07:00 "
                   "flow_cell_id=FBB01520";
  read_metadata_add_title(&m, t2, strlen(t2));
  test_eq_str(read_metadata_text(&m, READ_METADATA_RUN), "several");
  test_eq_str(read_metadata_text(&m, READ_METADATA_FLOW_CELL), "FBB01520");
  test_eq_str(read_metadata_text(&m, READ_METADATA_START),
              "2025-05-28T10:00:00+07:00");
  /* an Illumina read; a title which is neither */
  read_metadata_clear(&other);
  const char *t3 = "@A00123:45:HXYZ:1:1101:1000:2000 1:N:0:ACGT";
  read_metadata_add_title(&other, t3, strlen(t3));
  test_eq_str(read_metadata_text(&other, READ_METADATA_INSTRUMENT),
              "A00123");
  test_eq_str(read_metadata_text(&other, READ_METADATA_FLOW_CELL), "HXYZ");
  read_metadata_add_title(&other, "read1 length=100", 16);
  test_eq_str(read_metadata_text(&other, READ_METADATA_RUN), "45");
  /* merged, then saved and loaded */
  read_metadata_merge(&other, &m);
  test_eq_str(read_metadata_text(&other, READ_METADATA_FLOW_CELL),
              "several");
  char buf[2048], line[256];
  FILE *f = fmemopen(buf, sizeof(buf), "w+");
  read_metadata_save(&other, f);
  rewind(f);
  read_metadata_clear(&m);
  while (fgets(line, sizeof(line), f) != NULL)
    test_assert(read_metadata_load_line(&m, line));
  fclose(f);
  test_assert(memcmp(&m, &other, sizeof(read_metadata_t)) == 0);
  test_assert(!read_metadata_load_line(&m, "meta 99 0 x\n"));
}
