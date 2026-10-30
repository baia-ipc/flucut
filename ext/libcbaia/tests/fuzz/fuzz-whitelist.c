/* (c) 2026 Giorgio Gonnella */

/*
   Fuzzing harness: barcode whitelist. The whitelist is read and some
   sequences are searched in it, at both ends, with distances 0 to 2; the
   barcodes found exactly are counted and the cells identified.
*/

#include "fuzz-common.h"
#include "cell-caller.h"
#include "whitelist.h"

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  static bool initialized = false;
  static const char *queries[] = {"ACGTACGTACGTACGTACGT", "", "A",
    "NNNNNNNNNNNNNNNNNN", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"};
  if (!initialized) {
    fuzz_record_fds();
    initialized = true;
  }
  const char *path = fuzz_write_input(data, size);
  FUZZ_RUN({
    whitelist_t *wl = whitelist_new(path, NULL);
    cell_caller_t *cc = cell_caller_new();
    char out[WHITELIST_MAX_BARCODE_LENGTH + 1];
    for (size_t q = 0; q < sizeof(queries) / sizeof(queries[0]); q++)
      for (long d = 0; d <= 2; d++)
        for (int from_end = 0; from_end < 2; from_end++) {
          whitelist_match_t m;
          whitelist_match(wl, queries[q], strlen(queries[q]), from_end, d,
                          &m);
          if (m.status == WHITELIST_EXACT || m.status == WHITELIST_CORRECTED)
            whitelist_decode(wl, m.barcode, out);
          if (m.status == WHITELIST_EXACT)
            cell_caller_add(cc, m.barcode);
        }
    whitelist_t *cells = cell_caller_run(cc, 2,
        whitelist_barcode_length(wl), NULL, NULL, NULL);
    whitelist_delete(cells);
    cell_caller_delete(cc);
    whitelist_delete(wl);
  });
  return 0;
}
