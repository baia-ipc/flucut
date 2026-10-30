/* (c) 2026 Giorgio Gonnella */

/*
   Fuzzing harness: files of reads (FASTQ, optionally gzipped, or BAM,
   recognized by the reads reader). The records are read in small batches;
   each record is written again as FASTQ and as an unaligned BAM record,
   and its sequence is reverse complemented; the reads are also counted,
   and the completeness of the file is checked.
*/

#include "fuzz-common.h"
#include "fastq-writer.h"
#include "reads-reader.h"
#include "sequtil.h"
#include "strbuf.h"
#include "ubam-writer.h"

static fastq_batch_t batch;
static strbuf_t out;
static char *rc = NULL;
static size_t rc_size = 0;
static reads_reader_t *reader = NULL;

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  static bool initialized = false;
  if (!initialized) {
    fastq_batch_init(&batch);
    strbuf_init(&out);
    fuzz_record_fds();
    initialized = true;
  }
  const char *path = fuzz_write_input(data, size);
  (void)reads_reader_file_complete(path);
  FUZZ_RUN_CLEANUP({
    (void)reads_reader_count(path, NULL);
    reader = reads_reader_new(path, NULL);
    size_t n;
    while ((n = reads_reader_read_batch(reader, &batch, 3)) > 0)
      for (size_t i = 0; i < n; i++) {
        const fastq_record_t *r = batch.records + i;
        const char *d = batch.data.s;
        fastq_writer_record(&out, d + r->title_from, r->title_len,
                            d + r->seq_from, d + r->qual_from, r->seq_len);
        /* missing qualities (BAM input) written again as missing */
        ubam_writer_record(&out, d + r->title_from, r->title_len,
                           d + r->seq_from, r->qual_missing ? NULL :
                           d + r->qual_from, r->seq_len);
        if (r->seq_len + 1 > rc_size) {
          rc_size = r->seq_len + 1;
          rc = realloc(rc, rc_size);
          if (rc == NULL) {
            perror("realloc");
            abort();
          }
        }
        reverse_complement(d + r->seq_from, r->seq_len, rc);
        strbuf_clear(&out);
      }
    reads_reader_delete(reader);
    reader = NULL;
  }, {
    /* the reader buffers are large: free them also after an error */
    reads_reader_delete(reader);
    reader = NULL;
  });
  return 0;
}
