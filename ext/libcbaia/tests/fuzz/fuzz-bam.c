/* (c) 2026 Giorgio Gonnella */

/*
   Fuzzing harness: BAM files (the BGZF/gzip decompression is done by
   zlib). The header is encoded again as an unaligned BAM header, the
   records are read in small batches and encoded again as unaligned BAM
   records (with the tags of the input, and the missing qualities written
   again as missing); the records are also counted, and the end-of-file
   marker is checked.
*/

#include "fuzz-common.h"
#include "bam-reader.h"
#include "strbuf.h"
#include "ubam-writer.h"

static fastq_batch_t batch;
static strbuf_t out;
static bam_reader_t *reader = NULL;

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  static bool initialized = false;
  if (!initialized) {
    fastq_batch_init(&batch);
    strbuf_init(&out);
    fuzz_record_fds();
    initialized = true;
  }
  const char *path = fuzz_write_input(data, size);
  if (!bam_reader_is_bam(path, NULL))
    return 0;
  (void)bam_reader_has_eof_marker(path);
  FUZZ_RUN_CLEANUP({
    (void)bam_reader_count(path, NULL);
    reader = bam_reader_new(path, NULL);
    ubam_writer_header(&out, bam_reader_header(reader));
    size_t n;
    while ((n = bam_reader_read_batch(reader, &batch, 3)) > 0)
      for (size_t i = 0; i < n; i++) {
        const fastq_record_t *r = batch.records + i;
        const char *d = batch.data.s;
        /* missing qualities written again as missing (NULL) */
        ubam_writer_record(&out, d + r->title_from, r->title_len,
                           d + r->seq_from, r->qual_missing ? NULL :
                           d + r->qual_from, r->seq_len);
        strbuf_clear(&out);
      }
    bam_reader_delete(reader);
    reader = NULL;
  }, {
    bam_reader_delete(reader);
    reader = NULL;
  });
  return 0;
}
