/* (c) 2026 Giorgio Gonnella */

/*
   The platform of the reads (-p), read from the first reads of the first
   input file (-p auto, the default) or given, and in both cases checked
   against those reads before the trimming (in the live mode, at the
   first file ready), so that the settings of one platform are not used
   on the reads of the other without saying so.

   The two platforms give reads of different lengths, for physical
   reasons rather than by convention. An Illumina instrument reads a
   fixed number of cycles: no read is longer than that (151 or 251
   bases, at most 301), and most reads have exactly that length, the
   others being shorter where the instrument or a later step cut an
   adapter. Even two reads of a pair joined end to end are shorter than
   PLATFORM_MERGED_LONGEST. A nanopore read is a whole molecule: in a run
   of the RT-PCR most are whole amplicons of 900 to 2,400 bases, and
   their lengths have no common maximum. The read titles and the file
   names are not used: SRA, ENA and many pipelines rewrite them.

   From the lengths of the first PLATFORM_CHECK_READS reads:
   - illumina: none is longer than PLATFORM_ILLUMINA_LONGEST, and at
     least PLATFORM_FIXED_SHARE of them are within one base of the
     longest (reads of a fixed length);
   - nanopore: at least PLATFORM_LONG_SHARE of them are longer than
     PLATFORM_MERGED_LONGEST (longer than any Illumina read);
   - otherwise the evidence is weak: likely Illumina if none is longer
     than PLATFORM_MERGED_LONGEST, likely nanopore if some are;
   - too few reads to tell below PLATFORM_CHECK_MIN_READS.

   The thresholds lie in a wide gap between the two platforms, measured
   on 2026-11-28 on 86 Illumina files (43 paired runs of influenza A and
   B; MiSeq, NextSeq) and 344 nanopore files (MinION, GridION; R9.4.1
   and R10.4.1; 2022 to 2026): in each Illumina file the longest of the
   first 1,000 reads had 151 or 251 bases, and 51 to 100% of the reads
   were within one base of it; in each nanopore file the longest had at
   least 753 bases and at most 1.1% of the reads were within one base of
   it, and in each run of the RT-PCR at least 5% of the reads were
   longer than 600 bases. On blocks of 20, 100 and 1,000 consecutive
   reads (up to 20 per file), every Illumina block was called illumina,
   and no nanopore block was; at 1,000 reads 71 of 4,349 nanopore blocks
   were weak, all of four runs of random-primed cDNA (no RT-PCR, few
   reads longer than 600 bases). A first file of fewer reads gives more
   weak calls (on nanopore runs of the RT-PCR, none of the blocks of
   1,000 reads, 1.0% of those of 100, 1.5% of those of 20), which change
   no setting
   (../platform_check of the manuscript).

   What is done with it:
   - -p auto: the platform the reads clearly are; where the evidence is
     weak or there are too few reads, nanopore (the settings change only
     on clear evidence);
   - -p given: that platform. Where the reads clearly are of the other
     one, flucut stops and says what to pass instead, as the output
     would otherwise be quietly wrong (the nanopore settings cut
     fragments where a primer matches by chance inside them; the
     Illumina settings miss the primers behind the flanks of whole
     amplicons); --no-platform-check runs it anyway.
   Whenever the reads look like the other platform than the one used,
   clearly or weakly, a warning is written and the reports say so
   (the mismatch).
*/

#ifndef PLATFORM_CHECK_H
#define PLATFORM_CHECK_H

#include <stdbool.h>
#include <stddef.h>
#include "auto-settings.h"
#include "logger.h"
#include "options.h"
#include "trimmer.h"

#define PLATFORM_CHECK_READS 1000
#define PLATFORM_CHECK_MIN_READS 20
#define PLATFORM_ILLUMINA_LONGEST 320
#define PLATFORM_FIXED_SHARE 0.2
#define PLATFORM_MERGED_LONGEST 600
#define PLATFORM_LONG_SHARE 0.01

/* what the lengths of the reads say (see above) */
typedef enum {
  PLATFORM_READS_PENDING,          /* not read yet */
  PLATFORM_READS_FEW,              /* too few reads to tell */
  PLATFORM_READS_ILLUMINA,
  PLATFORM_READS_NANOPORE,
  PLATFORM_READS_LIKELY_ILLUMINA,
  PLATFORM_READS_LIKELY_NANOPORE,
} platform_reads_t;

/* the lengths of the reads looked at */
typedef struct {
  size_t n_reads;
  size_t longest;
  size_t at_longest;   /* reads within one base of the longest */
  size_t long_reads;   /* reads longer than PLATFORM_MERGED_LONGEST */
} platform_lengths_t;

typedef struct platform_check_t {
  bool automatic;      /* -p auto */
  bool stop;           /* stop if the reads are clearly of the other
                          platform (not --no-platform-check) */
  bool pending;
  bool illumina;       /* the platform used */
  platform_lengths_t lengths;
  platform_reads_t reads;
  /* the reads look like the other platform than the one used */
  bool mismatch;
  char *file;          /* the file looked at (owned) */
} platform_check_t;

/* Initialize the check from the options <opts> */
void platform_check_init(platform_check_t *self, const options_t *opts);

/* The lengths of the <n> reads of <lengths> */
platform_lengths_t platform_lengths_of(const size_t *lengths, size_t n);

/* What the lengths <l> say */
platform_reads_t platform_reads_classify(const platform_lengths_t *l);

/*
   The platform to use (illumina: true) for the reads <reads>, with -p
   auto if <automatic>, otherwise with the platform <given_illumina>;
   <mismatch> is set if the reads look like the other platform. Returns
   whether the reads clearly contradict the platform given.
*/
bool platform_check_choose(bool automatic, bool given_illumina,
    platform_reads_t reads, bool *illumina, bool *mismatch);

/*
   Decide the platform from the first reads of <file>, if still to be
   decided, and set its settings in <opts>, <config> and <autos> (the
   decision is reported to <logger>). Exits the program if the reads
   clearly contradict the platform given (unless --no-platform-check).
*/
void platform_check_decide(platform_check_t *self, options_t *opts,
    trimmer_config_t *config, auto_settings_t *autos, const char *file,
    Logger *logger);

/* The name of <reads> in the JSON summary (NULL: pending) */
const char* platform_reads_name(platform_reads_t reads);

/* What the lengths of the reads looked at say, as a sentence without
   its full stop, into <buf> of size <size> */
void platform_check_evidence(const platform_check_t *self, char *buf,
    size_t size);

/* Free the memory */
void platform_check_free(platform_check_t *self);

/* Unit tests */
void platform_check_test(void);

#endif
