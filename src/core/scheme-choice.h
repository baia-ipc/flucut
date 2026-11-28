/* (c) 2026 Giorgio Gonnella */

/*
   The choice of the primer scheme of influenza A from the reads (-s
   auto).

   The universal primers of influenza A all end in the same conserved
   termini and differ in the 5' tail in front of them (see primers.h):
   the protocols are told apart by their tails, not by the termini. A
   scheme which does not carry the tail of the protocol still finds an
   untailed terminus, but a tailed primer of another protocol is not
   found at all: its tail bases are all errors.

   The scheme is therefore chosen by the tails found in the first reads
   of the first input file, at each end of the amplicons apart (a
   protocol can use a tailed primer at one end and an untailed one at the
   other, as the default scheme uni12-mbt13 does):

   - at the start of each read and of its reverse complement (within
     the window of the primers and the length of a primer), the conserved
     terminus of each end is searched with at most
     SCHEME_CHOICE_TERMINUS_ERRORS errors (the terminus of influenza B,
     which differs from that of influenza A at two or three bases, is
     then not taken for it);
   - where it is found, each tail of the built-in schemes of influenza A
     for that end is aligned to the bases immediately in front of it:
     the tail must end where the terminus begins, with at most
     SCHEME_CHOICE_TAIL_ERRORS of its length in errors. A tail found
     anywhere else near the end of a read is not evidence of a primer
     (e.g. runs of G at the ends of the reads of two-colour Illumina
     chemistry, where no signal is read as G).

   The share of a tail is that of the read ends with the terminus which
   carry it: on Illumina reads, only a fragment beginning at an amplicon
   end carries a terminus, so a share over all the reads would mean
   nothing.

   An end is then
   - tailed with the leading tail, if it has at least SCHEME_CHOICE_MIN_SHARE
     of the read ends with the terminus and at least SCHEME_CHOICE_MIN_RATIO
     times the read ends of the next tail;
   - ambiguous, if it has the share but not the lead (two tails: a run
     holding more than one protocol);
   - untailed, if no tail has the share;
   - not decided, if fewer than SCHEME_CHOICE_MIN_ENDS read ends have
     the terminus.

   The scheme chosen is the built-in scheme of influenza A whose tails
   are those of the ends decided. The outcomes:
   - a scheme is chosen (untailed at both ends: uni);
   - ambiguous: an end has two tails; no single scheme fits;
   - not found: neither end has enough read ends with the terminus,
     which on these libraries usually means that the reads are not of
     this family of protocols;
   - no scheme: the tails found are not those of any built-in scheme.
   Where no scheme is chosen, the default scheme of influenza A
   (SCHEME_DEFAULT_A) is used, and the outcome is reported.

   The first SCHEME_CHOICE_READS reads are counted; unless both ends are
   then tailed, the reads up to SCHEME_CHOICE_MORE_READS are counted as
   well before deciding.

   The choice is within influenza A: the scheme of influenza B
   (SCHEME_DEFAULT_B) is the only one, and is added to the scheme chosen
   when the run may hold influenza B (--virus both); the primers of the
   type of a sample are then searched in its reads once its references
   say which type it is, as with any scheme of both types.
*/

#ifndef SCHEME_CHOICE_H
#define SCHEME_CHOICE_H

#include <stdbool.h>
#include <stddef.h>
#include "logger.h"
#include "options.h"
#include "primer-search.h"
#include "scheme.h"

#define SCHEME_CHOICE_READS 2000
#define SCHEME_CHOICE_MORE_READS 20000
#define SCHEME_CHOICE_TERMINUS_ERRORS 1
#define SCHEME_CHOICE_TAIL_ERRORS 0.2
#define SCHEME_CHOICE_MIN_ENDS 50
/* (a tail is not found in front of more than 0.6% of the read ends with
   the terminus in public runs without it, and in front of 81-98% at an
   end of the runs with it; in runs which mix a tailed and an untailed
   primer at an end, its share is anywhere from 13 to 66%: a lower
   share keeps those runs together, as having the tail) */
#define SCHEME_CHOICE_MIN_SHARE 0.05
#define SCHEME_CHOICE_MIN_RATIO 10
/* the most tails at an end, and the most schemes, of influenza A */
#define SCHEME_CHOICE_MAX_TAILS 8
#define SCHEME_CHOICE_MAX_SCHEMES 16
/* no tail (an untailed end of a scheme) */
#define SCHEME_CHOICE_NO_TAIL ((size_t)-1)

/* what the reads say about an end */
typedef enum {
  SCHEME_END_FEW,        /* too few read ends with the terminus */
  SCHEME_END_UNTAILED,   /* no tail has the share */
  SCHEME_END_TAILED,     /* the leading tail */
  SCHEME_END_AMBIGUOUS,  /* two tails */
} scheme_end_state_t;

/* the tails of an end and their counts */
typedef struct {
  size_t n_tails;
  const char *tail[SCHEME_CHOICE_MAX_TAILS];   /* (not 0-terminated) */
  size_t tail_len[SCHEME_CHOICE_MAX_TAILS];
  /* the read ends with the terminus, and those with each tail in front
     of it */
  size_t termini;
  size_t count[SCHEME_CHOICE_MAX_TAILS];
  scheme_end_state_t state;
  /* the leading tail and the next one (SCHEME_CHOICE_NO_TAIL: none) */
  size_t lead, second;
} scheme_end_counts_t;

typedef enum {
  SCHEME_CHOICE_CHOSEN,
  SCHEME_CHOICE_AMBIGUOUS,
  SCHEME_CHOICE_NOT_FOUND,
  SCHEME_CHOICE_NO_SCHEME,
} scheme_choice_outcome_t;

typedef struct {
  /* -s auto, and not yet decided */
  bool automatic, pending;
  /* the built-in schemes of influenza A and the tail of each at each
     end (an index in <ends>, SCHEME_CHOICE_NO_TAIL if untailed) */
  size_t n_schemes;
  const char *scheme[SCHEME_CHOICE_MAX_SCHEMES];
  size_t scheme_tail[SCHEME_CHOICE_MAX_SCHEMES][2];
  /* the left end (0) and the right end (1) */
  scheme_end_counts_t ends[2];
  scheme_choice_outcome_t outcome;
  /* the scheme of influenza A used: the one chosen, or the default */
  const char *chosen;
  size_t n_reads;
  char *file;   /* (owned) */
} scheme_choice_t;

/* Initialize from the options <opts>: automatic with -s auto, unless
   the run holds influenza B alone */
void scheme_choice_init(scheme_choice_t *self, const options_t *opts);

/* The scheme to load before the choice (the default of --virus with
   -s auto, the one of -s otherwise) */
const char* scheme_choice_initial(const options_t *opts);

/* Count the tails at the start of the read <seq> (length <len>) and at
   the start of its reverse complement, within <window> bases (see
   amplicon.h); <rc> has room for <len> bases */
void scheme_choice_count_read(scheme_choice_t *self,
    primer_search_t *search, const char *seq, size_t len,
    long window, char *rc);

/* Decide the state of each end and the outcome from the counts */
void scheme_choice_judge(scheme_choice_t *self);

/*
   If pending, choose the scheme from the first reads of <file>, and
   replace <scheme> by it (with the scheme of influenza B added if
   <virus> includes it); the choice is reported to <logger>
*/
void scheme_choice_decide(scheme_choice_t *self, scheme_t *scheme,
    virus_set_t virus, long window, const char *file, Logger *logger);

/* The name of an outcome, and of the state of an end (for the JSON
   summary) */
const char* scheme_choice_outcome_name(scheme_choice_outcome_t outcome);
const char* scheme_choice_end_state_name(scheme_end_state_t state);

/* What an outcome means, as text (for the log and the HTML report) */
const char* scheme_choice_outcome_text(scheme_choice_outcome_t outcome);

/* What the reads say about the end <e> (0: left, 1: right), as text */
void scheme_choice_end_text(const scheme_choice_t *self, size_t e,
    char *buf, size_t size);

/* Free the memory */
void scheme_choice_free(scheme_choice_t *self);

/* Unit tests */
void scheme_choice_test(void);

#endif
