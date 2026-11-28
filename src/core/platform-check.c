/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "alloc.h"
#include "fastq-reader.h"
#include "numutils.h"
#include "platform-check.h"
#include "reads-reader.h"
#include "strutils.h"
#include "test.h"

void platform_check_init(platform_check_t *self, const options_t *opts) {
  assert(self != NULL);
  assert(opts != NULL);
  *self = (platform_check_t){.automatic = opts->platform_auto,
      .stop = opts->platform_check, .pending = true,
      .illumina = opts->illumina, .reads = PLATFORM_READS_PENDING,
      .mismatch = false, .file = NULL};
}

platform_lengths_t platform_lengths_of(const size_t *lengths, size_t n) {
  platform_lengths_t l = {.n_reads = n};
  for (size_t i = 0; i < n; i++) {
    if (lengths[i] > l.longest)
      l.longest = lengths[i];
    if (lengths[i] > PLATFORM_MERGED_LONGEST)
      l.long_reads++;
  }
  for (size_t i = 0; i < n; i++)
    if (lengths[i] + 1 >= l.longest)
      l.at_longest++;
  return l;
}

platform_reads_t platform_reads_classify(const platform_lengths_t *l) {
  assert(l != NULL);
  double n = (double)l->n_reads;
  if (l->n_reads < PLATFORM_CHECK_MIN_READS)
    return PLATFORM_READS_FEW;
  if (l->longest <= PLATFORM_ILLUMINA_LONGEST &&
      (double)l->at_longest >= PLATFORM_FIXED_SHARE * n)
    return PLATFORM_READS_ILLUMINA;
  if (l->long_reads > 0 && (double)l->long_reads >= PLATFORM_LONG_SHARE * n)
    return PLATFORM_READS_NANOPORE;
  return l->long_reads == 0 ? PLATFORM_READS_LIKELY_ILLUMINA :
         PLATFORM_READS_LIKELY_NANOPORE;
}

bool platform_check_choose(bool automatic, bool given_illumina,
    platform_reads_t reads, bool *illumina, bool *mismatch) {
  bool leans_illumina = reads == PLATFORM_READS_ILLUMINA ||
                        reads == PLATFORM_READS_LIKELY_ILLUMINA;
  bool leans_nanopore = reads == PLATFORM_READS_NANOPORE ||
                        reads == PLATFORM_READS_LIKELY_NANOPORE;
  bool clear = reads == PLATFORM_READS_ILLUMINA ||
               reads == PLATFORM_READS_NANOPORE;
  /* -p auto: Illumina only on clear evidence */
  *illumina = automatic ? reads == PLATFORM_READS_ILLUMINA : given_illumina;
  *mismatch = (*illumina && leans_nanopore) ||
              (!*illumina && leans_illumina);
  return !automatic && *mismatch && clear;
}

/* the lengths of the first <max_reads> reads of <file> */
static platform_lengths_t measure(const char *file, size_t max_reads,
    Logger *logger) {
  size_t *lengths;
  XMALLOC(lengths, sizeof(size_t) * (max_reads > 0 ? max_reads : 1));
  reads_reader_t *reader = reads_reader_new(file, logger);
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  size_t n = 0, got;
  while (n < max_reads &&
         (got = reads_reader_read_batch(reader, &batch, max_reads - n)) > 0) {
    for (size_t i = 0; i < got; i++)
      lengths[n + i] = batch.records[i].seq_len;
    n += got;
  }
  fastq_batch_free(&batch);
  reads_reader_delete(reader);
  platform_lengths_t l = platform_lengths_of(lengths, n);
  free(lengths);
  return l;
}

/* the value of -p for a platform, and its name in a sentence */
static const char* platform_name(bool illumina) {
  return illumina ? "illumina" : "nanopore";
}

static const char* platform_text(bool illumina) {
  return illumina ? "Illumina" : "nanopore";
}

/* what the settings of a platform do to the reads of the other one */
static const char* consequence(bool illumina) {
  return illumina ?
      "the primers behind the flanks and adapters of whole amplicons would "
      "be missed and left in the reads" :
      "the reads would be cut where a primer matches by chance inside a "
      "fragment, removing template, and short fragments dropped";
}

void platform_check_evidence(const platform_check_t *self, char *buf,
    size_t size) {
  assert(self != NULL);
  const platform_lengths_t *l = &self->lengths;
  size_t n = l->n_reads;
  const char *file = self->file != NULL ? self->file : "";
  switch (self->reads) {
    case PLATFORM_READS_PENDING:
      snprintf(buf, size, "the reads are not read yet");
      break;
    case PLATFORM_READS_FEW:
      snprintf(buf, size, "%s has %zu reads, too few to tell the platform "
               "(at least %d needed)", file, n, PLATFORM_CHECK_MIN_READS);
      break;
    case PLATFORM_READS_ILLUMINA:
      snprintf(buf, size, "none of the first %zu reads of %s is longer than "
               "%zu bases, and %s are within one base of the longest: reads "
               "of a fixed length, as an Illumina run gives", n, file,
               l->longest, percent_text(l->at_longest, n).s);
      break;
    case PLATFORM_READS_NANOPORE:
      snprintf(buf, size, "%s of the first %zu reads of %s are longer than "
               "%d bases, which no Illumina read is: whole molecules, as a "
               "nanopore run gives", percent_text(l->long_reads, n).s, n,
               file, PLATFORM_MERGED_LONGEST);
      break;
    case PLATFORM_READS_LIKELY_ILLUMINA:
      snprintf(buf, size, "none of the first %zu reads of %s is longer than "
               "%d bases, as Illumina reads, but they are not of a fixed "
               "length (the longest has %zu bases, %s are within one base "
               "of it)", n, file, PLATFORM_MERGED_LONGEST, l->longest,
               percent_text(l->at_longest, n).s);
      break;
    case PLATFORM_READS_LIKELY_NANOPORE:
      snprintf(buf, size, "%zu of the first %zu reads of %s are longer than "
               "%d bases, which no Illumina read is, but they are fewer "
               "than %g%% of them", l->long_reads, n, file,
               PLATFORM_MERGED_LONGEST, 100 * PLATFORM_LONG_SHARE);
      break;
  }
}

/* set the settings of the platform chosen in <opts>, <config> and
   <autos> */
static void apply(bool illumina, options_t *opts, trimmer_config_t *config,
    auto_settings_t *autos) {
  if (illumina == opts->illumina)
    return;
  options_set_platform(opts, illumina);
  auto_settings_init(autos, opts);
  config->illumina = opts->illumina;
  config->window = opts->window;
  config->min_length = opts->min_length;
  config->adapters = opts->adapters;
  config->seeds_only = opts->illumina && !opts->exhaustive_search;
  config->max_errors = autos->max_errors.value;
}

/* report the decision */
static void log_decision(const platform_check_t *self, bool contradicts,
    Logger *logger) {
  char evidence[512];
  platform_check_evidence(self, evidence, sizeof(evidence));
  const char *used = platform_name(self->illumina);
  const char *other = platform_name(!self->illumina);
  const char *used_text = platform_text(self->illumina);
  const char *other_text = platform_text(!self->illumina);
  if (self->reads == PLATFORM_READS_FEW) {
    lprintf(logger, LOG_INFO, "Platform: %s, not checked: %s\n",
            used_text, evidence);
  } else if (!self->mismatch) {
    if (self->automatic)
      lprintf(logger, LOG_SUCCESS, "Platform: %s, read from the reads (-p "
              "auto): %s\n", used_text, evidence);
    else
      lprintf(logger, LOG_INFO, "Platform: %s (-p), as the reads look: "
              "%s\n", used_text, evidence);
  } else if (contradicts) {
    lprintf(logger, LOG_WARN, "the reads look like %s reads, not %s: %s. "
            "The %s settings are used as asked (-p %s, --no-platform-"
            "check), although with them %s\n", other_text, used_text,
            evidence, used_text, used, consequence(self->illumina));
  } else {
    lprintf(logger, LOG_WARN, "the reads may not be %s reads: %s. The %s "
            "settings are used%s; if the reads are %s reads, pass -p %s\n",
            used_text, evidence, used_text, self->automatic ?
            " (-p auto: Illumina only where the reads clearly are)" : "",
            other_text, other);
  }
}

void platform_check_decide(platform_check_t *self, options_t *opts,
    trimmer_config_t *config, auto_settings_t *autos, const char *file,
    Logger *logger) {
  assert(self != NULL);
  assert(opts != NULL);
  assert(config != NULL);
  assert(autos != NULL);
  assert(file != NULL);
  if (!self->pending)
    return;
  self->lengths = measure(file, PLATFORM_CHECK_READS, logger);
  self->reads = platform_reads_classify(&self->lengths);
  free(self->file);
  self->file = xstrdup(file);
  bool illumina;
  bool contradicts = platform_check_choose(self->automatic, opts->illumina,
      self->reads, &illumina, &self->mismatch);
  if (contradicts && self->stop) {
    char evidence[512];
    platform_check_evidence(self, evidence, sizeof(evidence));
    const char *given = platform_name(opts->illumina);
    LCHECK(false, logger, "the reads look like %s reads, not %s (-p %s): "
           "%s. With the %s settings, %s. Leave out -p, so that the "
           "platform is read from the reads, or pass -p %s; to use the %s "
           "settings anyway, add --no-platform-check.\n",
           platform_text(!opts->illumina), platform_text(opts->illumina),
           given, evidence, platform_text(opts->illumina),
           consequence(opts->illumina), platform_name(!opts->illumina),
           platform_text(opts->illumina));
  }
  apply(illumina, opts, config, autos);
  self->illumina = illumina;
  self->pending = false;
  log_decision(self, contradicts, logger);
}

const char* platform_reads_name(platform_reads_t reads) {
  switch (reads) {
    case PLATFORM_READS_PENDING: return NULL;
    case PLATFORM_READS_FEW: return "too_few_reads";
    case PLATFORM_READS_ILLUMINA: return "illumina";
    case PLATFORM_READS_NANOPORE: return "nanopore";
    case PLATFORM_READS_LIKELY_ILLUMINA: return "likely_illumina";
    case PLATFORM_READS_LIKELY_NANOPORE: return "likely_nanopore";
  }
  return NULL;
}

void platform_check_free(platform_check_t *self) {
  if (self == NULL)
    return;
  free(self->file);
  self->file = NULL;
}

/* the class of <n> reads, of which <n_fixed> have the length <fixed> and
   the others <other> */
static platform_reads_t classify_two(size_t n, size_t n_fixed,
    size_t fixed, size_t other) {
  size_t lengths[1000];
  assert(n <= 1000);
  for (size_t i = 0; i < n; i++)
    lengths[i] = i < n_fixed ? fixed : other;
  platform_lengths_t l = platform_lengths_of(lengths, n);
  return platform_reads_classify(&l);
}

void platform_check_test(void) {
  /* the lengths */
  {
    const size_t lengths[] = {151, 150, 149, 80, 151, 700, 1200};
    platform_lengths_t l = platform_lengths_of(lengths, 7);
    test_eq_size(l.n_reads, 7);
    test_eq_size(l.longest, 1200);
    test_eq_size(l.at_longest, 1);
    test_eq_size(l.long_reads, 2);
    l = platform_lengths_of(lengths, 5);
    test_eq_size(l.longest, 151);
    test_eq_size(l.at_longest, 3);   /* 151, 150, 151 */
    test_eq_size(l.long_reads, 0);
    l = platform_lengths_of(lengths, 0);
    test_eq_size(l.longest, 0);
    test_eq(platform_reads_classify(&l), PLATFORM_READS_FEW, int, "%d");
  }
  /* reads of a fixed length: Illumina, also with most of them shorter */
  test_eq(classify_two(1000, 1000, 151, 151), PLATFORM_READS_ILLUMINA,
          int, "%d");
  test_eq(classify_two(1000, 200, 301, 120), PLATFORM_READS_ILLUMINA,
          int, "%d");
  /* ... but not below the share, nor longer than an Illumina read */
  test_eq(classify_two(1000, 199, 301, 120),
          PLATFORM_READS_LIKELY_ILLUMINA, int, "%d");
  test_eq(classify_two(1000, 1000, PLATFORM_ILLUMINA_LONGEST + 1,
                       PLATFORM_ILLUMINA_LONGEST + 1),
          PLATFORM_READS_LIKELY_ILLUMINA, int, "%d");
  /* reads longer than an Illumina read: nanopore from 1% of them */
  test_eq(classify_two(1000, 10, 1500, 300), PLATFORM_READS_NANOPORE,
          int, "%d");
  test_eq(classify_two(1000, 9, 1500, 300),
          PLATFORM_READS_LIKELY_NANOPORE, int, "%d");
  test_eq(classify_two(1000, 9, PLATFORM_MERGED_LONGEST, 300),
          PLATFORM_READS_LIKELY_ILLUMINA, int, "%d");
  /* a few reads: one long read is enough, but not fewer than the
     minimum */
  test_eq(classify_two(PLATFORM_CHECK_MIN_READS, 1, 900, 400),
          PLATFORM_READS_NANOPORE, int, "%d");
  test_eq(classify_two(PLATFORM_CHECK_MIN_READS - 1, 19, 151, 151),
          PLATFORM_READS_FEW, int, "%d");
  /* the choice: -p auto */
  {
    bool illumina, mismatch;
    test_assert(!platform_check_choose(true, false, PLATFORM_READS_ILLUMINA,
                                       &illumina, &mismatch));
    test_assert(illumina && !mismatch);
    test_assert(!platform_check_choose(true, false, PLATFORM_READS_NANOPORE,
                                       &illumina, &mismatch));
    test_assert(!illumina && !mismatch);
    /* weak evidence: nanopore, with the mismatch where it leans to
       Illumina */
    test_assert(!platform_check_choose(true, false,
                                       PLATFORM_READS_LIKELY_ILLUMINA,
                                       &illumina, &mismatch));
    test_assert(!illumina && mismatch);
    test_assert(!platform_check_choose(true, false,
                                       PLATFORM_READS_LIKELY_NANOPORE,
                                       &illumina, &mismatch));
    test_assert(!illumina && !mismatch);
    test_assert(!platform_check_choose(true, false, PLATFORM_READS_FEW,
                                       &illumina, &mismatch));
    test_assert(!illumina && !mismatch);
  }
  /* the choice: a platform given */
  {
    bool illumina, mismatch;
    /* agreeing, or too few reads: as given */
    test_assert(!platform_check_choose(false, true, PLATFORM_READS_ILLUMINA,
                                       &illumina, &mismatch));
    test_assert(illumina && !mismatch);
    test_assert(!platform_check_choose(false, false, PLATFORM_READS_FEW,
                                       &illumina, &mismatch));
    test_assert(!illumina && !mismatch);
    test_assert(!platform_check_choose(false, true, PLATFORM_READS_FEW,
                                       &illumina, &mismatch));
    test_assert(illumina && !mismatch);
    /* clearly the other platform: a contradiction, both ways */
    test_assert(platform_check_choose(false, false, PLATFORM_READS_ILLUMINA,
                                      &illumina, &mismatch));
    test_assert(!illumina && mismatch);
    test_assert(platform_check_choose(false, true, PLATFORM_READS_NANOPORE,
                                      &illumina, &mismatch));
    test_assert(illumina && mismatch);
    /* weakly the other platform: the mismatch alone */
    test_assert(!platform_check_choose(false, false,
                                       PLATFORM_READS_LIKELY_ILLUMINA,
                                       &illumina, &mismatch));
    test_assert(!illumina && mismatch);
    test_assert(!platform_check_choose(false, true,
                                       PLATFORM_READS_LIKELY_NANOPORE,
                                       &illumina, &mismatch));
    test_assert(illumina && mismatch);
  }
  /* the settings of the platform chosen */
  {
    options_t opts = {.max_errors = OPTIONS_MAX_ERRORS_AUTO,
                      .platform_auto = true, .platform_check = true};
    options_set_platform(&opts, false);
    test_eq_long(opts.window, AMPLICON_DEFAULT_WINDOW);
    test_eq_size(opts.min_length, OPTIONS_DEFAULT_MIN_LENGTH);
    auto_settings_t autos;
    auto_settings_init(&autos, &opts);
    trimmer_config_t config = {.window = opts.window,
        .min_length = opts.min_length, .adapters = opts.adapters,
        .max_errors = autos.max_errors.value};
    apply(true, &opts, &config, &autos);
    test_assert(opts.illumina && config.illumina && config.seeds_only);
    test_eq_long(config.window, OPTIONS_ILLUMINA_WINDOW);
    test_eq_size(config.min_length, OPTIONS_ILLUMINA_MIN_LENGTH);
    test_eq_str(config.adapters, "none");
    test_eq(config.max_errors, AUTO_MAX_ERRORS_ILLUMINA, double, "%g");
    test_assert(!auto_settings_pending(&autos));
    auto_settings_free(&autos);
    /* the settings given are kept */
    options_t given = {.max_errors = 0.1, .window = 30,
                       .window_given = true, .platform_auto = true};
    options_set_platform(&given, false);
    test_eq_long(given.window, 30);
    auto_settings_init(&autos, &given);
    config = (trimmer_config_t){.window = given.window};
    apply(true, &given, &config, &autos);
    test_eq_long(config.window, 30);
    test_eq(config.max_errors, 0.1, double, "%g");
    auto_settings_free(&autos);
  }
}
