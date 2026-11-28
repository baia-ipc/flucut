/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "alloc.h"
#include "fastq-reader.h"
#include "iupac.h"
#include "primers.h"
#include "reads-reader.h"
#include "scheme-choice.h"
#include "sequtil.h"
#include "strutils.h"
#include "test.h"

/* the conserved terminus of each end of influenza A (the degenerate
   form of the left one covers both forms of its fourth base) */
static const char *const end_terminus[2] = {PRIMERS_TERMINUS_12_R,
                                            PRIMERS_TERMINUS_13};
static const char *const end_name[2] = {"left", "right"};

const char* scheme_choice_initial(const options_t *opts) {
  assert(opts != NULL);
  return strcmp(opts->scheme, OPTIONS_SCHEME_AUTO) == 0 ?
         scheme_default_for(opts->virus) : opts->scheme;
}

void scheme_choice_init(scheme_choice_t *self, const options_t *opts) {
  assert(self != NULL);
  assert(opts != NULL);
  memset(self, 0, sizeof(*self));
  self->automatic = strcmp(opts->scheme, OPTIONS_SCHEME_AUTO) == 0 &&
                    (opts->virus & VIRUS_SET_A) != 0;
  self->pending = self->automatic;
  self->chosen = SCHEME_DEFAULT_A;
  self->outcome = SCHEME_CHOICE_CHOSEN;
}

/* the tail of the built-in primer named <name>: the bases in front of
   its terminus, as it is ordered */
static const char* primer_tail(const char *name, size_t *len) {
  const primer_t *p = primers_builtin(name);
  assert(p != NULL);
  *len = p->terminus_from;
  return p->seq;
}

/* the index of the tail <tail> (length <len>) among those of the end
   <e>, added if new; SCHEME_CHOICE_NO_TAIL for no tail */
static size_t tail_index(scheme_choice_t *self, size_t e, const char *tail,
    size_t len) {
  if (len == 0)
    return SCHEME_CHOICE_NO_TAIL;
  scheme_end_counts_t *end = self->ends + e;
  for (size_t t = 0; t < end->n_tails; t++)
    if (end->tail_len[t] == len && memcmp(end->tail[t], tail, len) == 0)
      return t;
  assert(end->n_tails < SCHEME_CHOICE_MAX_TAILS);
  end->tail[end->n_tails] = tail;
  end->tail_len[end->n_tails] = len;
  return end->n_tails++;
}

/* the tail shared by the <n> primers <p> (as found in a read: those of
   the right end reverse-complemented; their tail is taken from the
   primer as it is ordered); false if they do not share one */
static bool shared_tail(const primer_t *p, size_t n, const char **tail,
    size_t *len) {
  if (n == 0)
    return false;
  *tail = primer_tail(p[0].name, len);
  for (size_t i = 1; i < n; i++) {
    size_t other_len;
    const char *other = primer_tail(p[i].name, &other_len);
    if (other_len != *len || memcmp(other, *tail, *len) != 0)
      return false;
  }
  return true;
}

/* the built-in schemes of influenza A whose primers share a tail at
   each end, and their tails */
static void collect_schemes(scheme_choice_t *self, Logger *logger) {
  if (self->n_schemes > 0)
    return;
  for (const char *const *n = scheme_builtin_names(); *n != NULL; n++) {
    scheme_t *s = scheme_new_builtin(*n, logger);
    const char *tail[2];
    size_t len[2];
    if (scheme_types(s) == VIRUS_SET_A &&
        shared_tail(s->left, s->n_left, tail + 0, len + 0) &&
        shared_tail(s->right, s->n_right, tail + 1, len + 1)) {
      assert(self->n_schemes < SCHEME_CHOICE_MAX_SCHEMES);
      size_t k = self->n_schemes++;
      self->scheme[k] = *n;
      for (size_t e = 0; e < 2; e++)
        self->scheme_tail[k][e] = tail_index(self, e, tail[e], len[e]);
    }
    scheme_delete(s);
  }
}

/* whether the characters <a> and <b> are equal (IUPAC codes) */
static bool chars_equal(const iupac_equal_chars_t *eq, char a, char b) {
  const iupac_equal_chars_t *e = eq + (unsigned char)a;
  for (unsigned char i = 0; i < e->n; i++)
    if (e->chars[i] == (unsigned char)b)
      return true;
  return false;
}

/*
   The errors of the tail <tail> (length <m>) aligned to the bases of
   <seq> which end at <to> (exclusive): the tail must end there, where
   the terminus begins; it may begin anywhere before (edit distance)
*/
#define TAIL_MAX_LEN 32
static long tail_errors(const char *tail, size_t m, const char *seq,
    size_t to) {
  assert(m <= TAIL_MAX_LEN);
  const iupac_equal_chars_t *eq = iupac_equal_chars();
  /* the read bases considered: as many as the tail and as many again
     (insertions beyond that cost more than the errors allowed) */
  size_t n = to < 2 * m ? to : 2 * m;
  /* row i: the last i bases of the tail; column j: the last j bases
     before <to> */
  long prev[2 * TAIL_MAX_LEN + 1], cur[2 * TAIL_MAX_LEN + 1];
  for (size_t j = 0; j <= n; j++)
    prev[j] = (long)j;
  for (size_t i = 1; i <= m; i++) {
    cur[0] = (long)i;
    char t = tail[m - i];
    for (size_t j = 1; j <= n; j++) {
      long best = prev[j - 1] + (chars_equal(eq, t, seq[to - j]) ? 0 : 1);
      if (prev[j] + 1 < best)
        best = prev[j] + 1;
      if (cur[j - 1] + 1 < best)
        best = cur[j - 1] + 1;
      cur[j] = best;
    }
    memcpy(prev, cur, sizeof(long) * (n + 1));
  }
  long best = prev[0];
  for (size_t j = 1; j <= n; j++)
    if (prev[j] < best)
      best = prev[j];
  return best;
}

/* count the termini and the tails at the start of <seq> */
static void count_start(scheme_choice_t *self, primer_search_t *search,
    const char *seq, size_t len, long window) {
  for (size_t e = 0; e < 2; e++) {
    scheme_end_counts_t *end = self->ends + e;
    size_t terminus_len = strlen(end_terminus[e]), longest = 0;
    for (size_t t = 0; t < end->n_tails; t++)
      if (end->tail_len[t] > longest)
        longest = end->tail_len[t];
    size_t region = (size_t)window + longest + terminus_len;
    primer_hit_t hit;
    if (!primer_search_run_seq(search, end_terminus[e], terminus_len, seq,
                               len < region ? len : region,
                               SCHEME_CHOICE_TERMINUS_ERRORS, false, false,
                               &hit))
      continue;
    end->termini++;
    size_t best = SCHEME_CHOICE_NO_TAIL;
    long best_errors = 0;
    for (size_t t = 0; t < end->n_tails; t++) {
      long errors = tail_errors(end->tail[t], end->tail_len[t], seq,
                                (size_t)hit.start);
      if (errors <= (long)floor(SCHEME_CHOICE_TAIL_ERRORS *
                                (double)end->tail_len[t]) &&
          (best == SCHEME_CHOICE_NO_TAIL || errors < best_errors)) {
        best = t;
        best_errors = errors;
      }
    }
    if (best != SCHEME_CHOICE_NO_TAIL)
      end->count[best]++;
  }
}

void scheme_choice_count_read(scheme_choice_t *self,
    primer_search_t *search, const char *seq, size_t len, long window,
    char *rc) {
  assert(self != NULL);
  count_start(self, search, seq, len, window);
  reverse_complement(seq, len, rc);
  count_start(self, search, rc, len, window);
}

/* the state of an end from its counts */
static void judge_end(scheme_end_counts_t *end) {
  end->lead = end->second = SCHEME_CHOICE_NO_TAIL;
  for (size_t t = 0; t < end->n_tails; t++) {
    if (end->lead == SCHEME_CHOICE_NO_TAIL ||
        end->count[t] > end->count[end->lead]) {
      end->second = end->lead;
      end->lead = t;
    } else if (end->second == SCHEME_CHOICE_NO_TAIL ||
               end->count[t] > end->count[end->second])
      end->second = t;
  }
  size_t lead = end->lead == SCHEME_CHOICE_NO_TAIL ? 0 :
                end->count[end->lead];
  size_t second = end->second == SCHEME_CHOICE_NO_TAIL ? 0 :
                  end->count[end->second];
  if (end->termini < SCHEME_CHOICE_MIN_ENDS)
    end->state = SCHEME_END_FEW;
  else if (lead == 0 ||
           (double)lead < SCHEME_CHOICE_MIN_SHARE * (double)end->termini)
    end->state = SCHEME_END_UNTAILED;
  else if ((double)lead < SCHEME_CHOICE_MIN_RATIO * (double)second)
    end->state = SCHEME_END_AMBIGUOUS;
  else
    end->state = SCHEME_END_TAILED;
}

/* whether the tail <tail> of a scheme at an end fits what the reads
   say there */
static bool end_fits(const scheme_end_counts_t *end, size_t tail) {
  switch (end->state) {
    case SCHEME_END_FEW:
      return true;
    case SCHEME_END_UNTAILED:
      return tail == SCHEME_CHOICE_NO_TAIL;
    case SCHEME_END_TAILED:
      return tail == end->lead;
    default:
      return false;
  }
}

void scheme_choice_judge(scheme_choice_t *self) {
  assert(self != NULL);
  for (size_t e = 0; e < 2; e++)
    judge_end(self->ends + e);
  self->chosen = SCHEME_DEFAULT_A;
  if (self->ends[0].state == SCHEME_END_FEW &&
      self->ends[1].state == SCHEME_END_FEW) {
    self->outcome = SCHEME_CHOICE_NOT_FOUND;
    return;
  }
  if (self->ends[0].state == SCHEME_END_AMBIGUOUS ||
      self->ends[1].state == SCHEME_END_AMBIGUOUS) {
    self->outcome = SCHEME_CHOICE_AMBIGUOUS;
    return;
  }
  /* a scheme is chosen by a tail: without one, the default is kept (see
     scheme-choice.h on why not the scheme without tails) */
  if (self->ends[0].state != SCHEME_END_TAILED &&
      self->ends[1].state != SCHEME_END_TAILED) {
    self->outcome = SCHEME_CHOICE_UNTAILED;
    return;
  }
  /* the schemes which fit both ends: one is chosen; of several (an end
     not decided), the default if it is one of them */
  size_t n_fit = 0, fit = 0;
  bool default_fits = false;
  for (size_t k = 0; k < self->n_schemes; k++) {
    if (!end_fits(self->ends + 0, self->scheme_tail[k][0]) ||
        !end_fits(self->ends + 1, self->scheme_tail[k][1]))
      continue;
    n_fit++;
    fit = k;
    if (strcmp(self->scheme[k], SCHEME_DEFAULT_A) == 0)
      default_fits = true;
  }
  if (n_fit == 1)
    self->chosen = self->scheme[fit];
  self->outcome = n_fit == 1 || default_fits ? SCHEME_CHOICE_CHOSEN :
                  SCHEME_CHOICE_NO_SCHEME;
}

const char* scheme_choice_outcome_name(scheme_choice_outcome_t outcome) {
  switch (outcome) {
    case SCHEME_CHOICE_CHOSEN:
      return "chosen";
    case SCHEME_CHOICE_UNTAILED:
      return "untailed";
    case SCHEME_CHOICE_AMBIGUOUS:
      return "ambiguous";
    case SCHEME_CHOICE_NOT_FOUND:
      return "not_found";
    default:
      return "no_scheme";
  }
}

const char* scheme_choice_end_state_name(scheme_end_state_t state) {
  switch (state) {
    case SCHEME_END_FEW:
      return "few";
    case SCHEME_END_UNTAILED:
      return "untailed";
    case SCHEME_END_TAILED:
      return "tailed";
    default:
      return "ambiguous";
  }
}

const char* scheme_choice_outcome_text(scheme_choice_outcome_t outcome) {
  switch (outcome) {
    case SCHEME_CHOICE_CHOSEN:
      return "chosen by the tails in front of the conserved termini";
    case SCHEME_CHOICE_UNTAILED:
      return "no tail in front of the conserved termini: the default is "
             "used";
    case SCHEME_CHOICE_AMBIGUOUS:
      return "two tails at an end, so no single scheme fits: the default "
             "is used";
    case SCHEME_CHOICE_NOT_FOUND:
      return "too few read ends with a conserved terminus: the default is "
             "used";
    default:
      return "no built-in scheme has the tails found: the default is used";
  }
}

void scheme_choice_end_text(const scheme_choice_t *self, size_t e,
    char *buf, size_t size) {
  assert(self != NULL);
  assert(e < 2);
  const scheme_end_counts_t *end = self->ends + e;
  size_t lead = end->lead, second = end->second;
  int lead_len = lead == SCHEME_CHOICE_NO_TAIL ? 0 :
                 (int)end->tail_len[lead];
  int second_len = second == SCHEME_CHOICE_NO_TAIL ? 0 :
                   (int)end->tail_len[second];
  const char *lead_tail = lead == SCHEME_CHOICE_NO_TAIL ? "" :
                          end->tail[lead];
  const char *second_tail = second == SCHEME_CHOICE_NO_TAIL ? "" :
                            end->tail[second];
  size_t lead_n = lead == SCHEME_CHOICE_NO_TAIL ? 0 : end->count[lead];
  size_t second_n = second == SCHEME_CHOICE_NO_TAIL ? 0 :
                    end->count[second];
  switch (end->state) {
    case SCHEME_END_FEW:
      snprintf(buf, size, "%s end: %zu read ends with the terminus, too "
               "few", end_name[e], end->termini);
      break;
    case SCHEME_END_UNTAILED:
      if (lead_n == 0)
        snprintf(buf, size, "%s end: no tail in %zu read ends with the "
                 "terminus", end_name[e], end->termini);
      else
        snprintf(buf, size, "%s end: no tail (the most frequent, %.*s, in "
                 "%zu of %zu read ends with the terminus)", end_name[e],
                 lead_len, lead_tail, lead_n, end->termini);
      break;
    case SCHEME_END_TAILED:
      snprintf(buf, size, "%s end: tail %.*s in %zu of %zu read ends with "
               "the terminus (next: %.*s, %zu)", end_name[e], lead_len,
               lead_tail, lead_n, end->termini, second_len, second_tail,
               second_n);
      break;
    default:
      snprintf(buf, size, "%s end: two tails, %.*s in %zu and %.*s in %zu "
               "of %zu read ends with the terminus", end_name[e], lead_len,
               lead_tail, lead_n, second_len, second_tail, second_n,
               end->termini);
  }
}

/* report the choice of <file> */
static void log_choice(const scheme_choice_t *self, const char *name,
    const char *file, Logger *logger) {
  char ends[2][256];
  for (size_t e = 0; e < 2; e++)
    scheme_choice_end_text(self, e, ends[e], sizeof(ends[e]));
  bool used = self->outcome == SCHEME_CHOICE_CHOSEN ||
              self->outcome == SCHEME_CHOICE_UNTAILED;
  lprintf(logger, used ? LOG_INFO : LOG_WARN, "Primer scheme %s (-s auto): "
          "%s, in the first %zu reads of %s; %s; %s\n", name,
          scheme_choice_outcome_text(self->outcome), self->n_reads, file,
          ends[0], ends[1]);
}

void scheme_choice_decide(scheme_choice_t *self, scheme_t *scheme,
    virus_set_t virus, long window, const char *file, Logger *logger) {
  assert(self != NULL);
  assert(scheme != NULL);
  assert(file != NULL);
  if (!self->pending)
    return;
  clock_t started = clock();
  collect_schemes(self, logger);
  reads_reader_t *reader = reads_reader_new(file, logger);
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  primer_search_t *search = primer_search_new();
  char *rc = NULL;
  size_t rc_size = 0, n = 0;
  const size_t stages[] = {SCHEME_CHOICE_READS, SCHEME_CHOICE_MORE_READS};
  for (size_t stage = 0; stage < 2; stage++) {
    size_t got;
    while (n < stages[stage] &&
           (got = reads_reader_read_batch(reader, &batch,
                                          stages[stage] - n)) > 0) {
      for (size_t i = 0; i < got; i++) {
        const fastq_record_t *rec = batch.records + i;
        if (rec->seq_len + 1 > rc_size) {
          rc_size = rec->seq_len + 1;
          XREALLOC(rc, rc_size);
        }
        scheme_choice_count_read(self, search, batch.data.s + rec->seq_from,
                                 rec->seq_len, window, rc);
      }
      n += got;
    }
    scheme_choice_judge(self);
    /* both ends tailed, or no more reads */
    if ((self->outcome == SCHEME_CHOICE_CHOSEN &&
         self->ends[0].state == SCHEME_END_TAILED &&
         self->ends[1].state == SCHEME_END_TAILED) || n < stages[stage])
      break;
  }
  free(rc);
  primer_search_delete(search);
  fastq_batch_free(&batch);
  reads_reader_delete(reader);
  self->n_reads = n;
  free(self->file);
  self->file = xstrdup(file);
  self->pending = false;
  /* the scheme of influenza B is added if the run may hold it */
  char name[128];
  snprintf(name, sizeof(name), "%s%s", self->chosen,
           (virus & VIRUS_SET_B) != 0 ? "+" SCHEME_DEFAULT_B : "");
  scheme_replace(scheme, scheme_new_builtin(name, logger));
  log_choice(self, scheme->name, file, logger);
  lprintf(logger, LOG_DEBUG, "Primer scheme chosen in %.3f s of processor "
          "time\n", (double)(clock() - started) / CLOCKS_PER_SEC);
}

void scheme_choice_free(scheme_choice_t *self) {
  if (self == NULL)
    return;
  free(self->file);
  self->file = NULL;
}

/* the index of the scheme <name> among those of <c> */
static size_t scheme_index(const scheme_choice_t *c, const char *name) {
  for (size_t k = 0; k < c->n_schemes; k++)
    if (strcmp(c->scheme[k], name) == 0)
      return k;
  return SIZE_MAX;
}

/* set the counts of the end <e> of <c>: <termini> and those of the tail
   of the scheme <a> and of the scheme <b> (NULL: none) */
static void set_counts(scheme_choice_t *c, size_t e, size_t termini,
    const char *a, size_t count_a, const char *b, size_t count_b) {
  scheme_end_counts_t *end = c->ends + e;
  end->termini = termini;
  memset(end->count, 0, sizeof(end->count));
  end->count[c->scheme_tail[scheme_index(c, a)][e]] = count_a;
  if (b != NULL)
    end->count[c->scheme_tail[scheme_index(c, b)][e]] = count_b;
}

/* a read of an amplicon of the primers <left> and <right> (as ordered):
   <flank>, the left primer, <template>, the reverse complement of the
   right primer, <flank> */
static char* test_read(const char *flank, const char *left,
    const char *template, const char *right) {
  size_t len = 2 * strlen(flank) + strlen(left) + strlen(template) +
               strlen(right);
  char *r;
  XMALLOC(r, len + 1);
  char *rc_right;
  XMALLOC(rc_right, strlen(right) + 1);
  reverse_complement(right, strlen(right), rc_right);
  rc_right[strlen(right)] = '\0';
  snprintf(r, len + 1, "%s%s%s%s%s", flank, left, template, rc_right,
           flank);
  free(rc_right);
  return r;
}

void scheme_choice_test(void) {
  Logger *logger = new_logger(stderr, LOG_ERROR);
  scheme_choice_t c;
  options_t opts;
  memset(&opts, 0, sizeof(opts));
  opts.scheme = OPTIONS_SCHEME_AUTO;
  opts.virus = VIRUS_SET_B;
  scheme_choice_init(&c, &opts);
  test_assert(!c.automatic);
  test_eq_str(scheme_choice_initial(&opts), SCHEME_DEFAULT_B);
  opts.virus = VIRUS_SET_BOTH;
  scheme_choice_init(&c, &opts);
  test_assert(c.automatic && c.pending);
  test_eq_str(scheme_choice_initial(&opts), SCHEME_DEFAULT_BOTH);
  opts.scheme = "opti";
  test_eq_str(scheme_choice_initial(&opts), "opti");
  /* the schemes of influenza A and their tails */
  collect_schemes(&c, logger);
  size_t def = scheme_index(&c, SCHEME_DEFAULT_A),
         mbt = scheme_index(&c, "mbt"), uni = scheme_index(&c, "uni"),
         opti = scheme_index(&c, "opti"), inf = scheme_index(&c, "uni-inf");
  test_assert(def != SIZE_MAX && mbt != SIZE_MAX && uni != SIZE_MAX &&
              opti != SIZE_MAX && inf != SIZE_MAX);
  test_eq_size(scheme_index(&c, SCHEME_DEFAULT_B), SIZE_MAX);
  test_eq_size(c.scheme_tail[def][0], SCHEME_CHOICE_NO_TAIL);
  test_eq_size(c.scheme_tail[uni][0], SCHEME_CHOICE_NO_TAIL);
  test_eq_size(c.scheme_tail[uni][1], SCHEME_CHOICE_NO_TAIL);
  test_eq_size(c.scheme_tail[def][1], c.scheme_tail[mbt][1]);
  size_t t = c.scheme_tail[opti][0];
  test_eq_mem(c.ends[0].tail[t], c.ends[0].tail_len[t], PRIMERS_OPTI_TAIL);
  t = c.scheme_tail[inf][1];
  test_eq_mem(c.ends[1].tail[t], c.ends[1].tail_len[t],
              PRIMERS_INF_TAIL_RIGHT);
  /* the mbt tail overlaps the terminus by two bases: the tail is what
     lies in front of the terminus */
  t = c.scheme_tail[mbt][0];
  test_eq_mem(c.ends[0].tail[t], c.ends[0].tail_len[t], "ACGCGTGATC");
  /* the decision: a tail at both ends */
  set_counts(&c, 0, 1990, "opti", 1930, "uni-inf", 2);
  set_counts(&c, 1, 1980, "opti", 1900, NULL, 0);
  scheme_choice_judge(&c);
  test_eq_int((int)c.outcome, SCHEME_CHOICE_CHOSEN);
  test_eq_str(c.chosen, "opti");
  /* a short tail, often lost to errors, still has the share */
  set_counts(&c, 0, 1900, "uni-inf", 1200, "mbt", 1);
  set_counts(&c, 1, 1900, "uni-inf", 1500, NULL, 0);
  scheme_choice_judge(&c);
  test_eq_str(c.chosen, "uni-inf");
  /* untailed at the left end, the tail of mbt at the right end: the
     default */
  set_counts(&c, 0, 1900, "uni-inf", 3, NULL, 0);
  set_counts(&c, 1, 1900, "mbt", 1800, NULL, 0);
  scheme_choice_judge(&c);
  test_eq_int((int)c.outcome, SCHEME_CHOICE_CHOSEN);
  test_eq_int((int)c.ends[0].state, SCHEME_END_UNTAILED);
  test_eq_str(c.chosen, SCHEME_DEFAULT_A);
  /* untailed at both ends: the default, not the scheme without tails */
  set_counts(&c, 1, 1900, "mbt", 5, NULL, 0);
  scheme_choice_judge(&c);
  test_eq_int((int)c.outcome, SCHEME_CHOICE_UNTAILED);
  test_eq_str(c.chosen, SCHEME_DEFAULT_A);
  /* ... also with too few at the other end */
  set_counts(&c, 0, 30, "opti", 1, NULL, 0);
  scheme_choice_judge(&c);
  test_eq_int((int)c.outcome, SCHEME_CHOICE_UNTAILED);
  test_eq_str(c.chosen, SCHEME_DEFAULT_A);
  /* two tails at an end: ambiguous, the default */
  set_counts(&c, 0, 1900, "opti", 900, "mbt", 850);
  set_counts(&c, 1, 1900, "opti", 950, "mbt", 800);
  scheme_choice_judge(&c);
  test_eq_int((int)c.outcome, SCHEME_CHOICE_AMBIGUOUS);
  test_eq_int((int)c.ends[0].state, SCHEME_END_AMBIGUOUS);
  test_eq_str(c.chosen, SCHEME_DEFAULT_A);
  /* ... also when the second has much less, but not a tenth */
  set_counts(&c, 0, 1900, "opti", 1500, "mbt", 160);
  set_counts(&c, 1, 1900, "opti", 1500, NULL, 0);
  scheme_choice_judge(&c);
  test_eq_int((int)c.outcome, SCHEME_CHOICE_AMBIGUOUS);
  /* too few termini at both ends: not found */
  set_counts(&c, 0, SCHEME_CHOICE_MIN_ENDS - 1, "opti", 40, NULL, 0);
  set_counts(&c, 1, 3, "opti", 3, NULL, 0);
  scheme_choice_judge(&c);
  test_eq_int((int)c.outcome, SCHEME_CHOICE_NOT_FOUND);
  test_eq_str(c.chosen, SCHEME_DEFAULT_A);
  /* too few at one end: the other one decides, if it tells one scheme */
  set_counts(&c, 1, 400, "opti", 390, NULL, 0);
  scheme_choice_judge(&c);
  test_eq_str(c.chosen, "opti");
  set_counts(&c, 1, 400, "mbt", 390, NULL, 0);
  scheme_choice_judge(&c);
  test_eq_int((int)c.outcome, SCHEME_CHOICE_CHOSEN);
  test_eq_str(c.chosen, SCHEME_DEFAULT_A);
  /* tails of no built-in scheme */
  set_counts(&c, 0, 1900, "opti", 1800, NULL, 0);
  set_counts(&c, 1, 1900, "mbt", 1800, NULL, 0);
  scheme_choice_judge(&c);
  test_eq_int((int)c.outcome, SCHEME_CHOICE_NO_SCHEME);
  test_eq_str(c.chosen, SCHEME_DEFAULT_A);
  /* the counting: a read of the Opti1 primers, in either orientation */
  for (size_t e = 0; e < 2; e++) {
    c.ends[e].termini = 0;
    memset(c.ends[e].count, 0, sizeof(c.ends[e].count));
  }
  primer_search_t *search = primer_search_new();
  const char *template = "TTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAA";
  char *r = test_read("ACGTTGCAACGTTGCA", PRIMERS_OPTI_TAIL
                      PRIMERS_TERMINUS_12, template,
                      PRIMERS_OPTI_TAIL PRIMERS_TERMINUS_13);
  char rc[256];
  scheme_choice_count_read(&c, search, r, strlen(r), 100, rc);
  free(r);
  test_eq_size(c.ends[0].termini, 1);
  test_eq_size(c.ends[1].termini, 1);
  test_eq_size(c.ends[0].count[c.scheme_tail[opti][0]], 1);
  test_eq_size(c.ends[1].count[c.scheme_tail[opti][1]], 1);
  /* a tail with an error in front of the terminus counts; one not
     directly in front of it (three bases between) does not, at either
     end */
  r = test_read("ACGTTGCAACGTTGCA", "GTTACGAGCC" PRIMERS_TERMINUS_12,
                template, PRIMERS_INF_TAIL_RIGHT "TCA" PRIMERS_TERMINUS_13);
  scheme_choice_count_read(&c, search, r, strlen(r), 100, rc);
  free(r);
  test_eq_size(c.ends[0].count[c.scheme_tail[opti][0]], 2);
  test_eq_size(c.ends[1].termini, 2);
  for (size_t k = 0; k < c.ends[1].n_tails; k++)
    test_eq_size(c.ends[1].count[k], k == c.scheme_tail[opti][1] ? 1 : 0);
  r = test_read("ACGTTGCAACGTTGCA", "GGGGGGTCA" PRIMERS_TERMINUS_12,
                template, PRIMERS_TERMINUS_13);
  scheme_choice_count_read(&c, search, r, strlen(r), 100, rc);
  free(r);
  test_eq_size(c.ends[0].termini, 3);
  test_eq_size(c.ends[0].count[c.scheme_tail[inf][0]], 0);
  test_eq_size(c.ends[1].termini, 3);
  /* an untailed terminus, and the termini of influenza B, which are not
     taken for those of influenza A */
  r = test_read("ACGTTGCAACGTTGCA", PRIMERS_TERMINUS_12, template,
                PRIMERS_IBV_TAIL_RIGHT PRIMERS_TERMINUS_B_10 "CGAGC");
  scheme_choice_count_read(&c, search, r, strlen(r), 100, rc);
  free(r);
  test_eq_size(c.ends[0].termini, 4);
  test_eq_size(c.ends[0].count[c.scheme_tail[opti][0]], 2);
  test_eq_size(c.ends[1].termini, 3);
  r = test_read("ACGTTGCAACGTTGCA",
                PRIMERS_IBV_TAIL_LEFT PRIMERS_TERMINUS_B_9 "GGAGC", template,
                PRIMERS_TERMINUS_13);
  scheme_choice_count_read(&c, search, r, strlen(r), 100, rc);
  free(r);
  test_eq_size(c.ends[0].termini, 4);
  test_eq_size(c.ends[1].termini, 4);
  /* a terminus beyond the window is not counted */
  r = test_read("ACGTTGCAACGTTGCA", PRIMERS_TERMINUS_12, template,
                PRIMERS_TERMINUS_13);
  scheme_choice_count_read(&c, search, r, strlen(r), 0, rc);
  free(r);
  test_eq_size(c.ends[0].termini, 4);
  test_eq_size(c.ends[1].termini, 4);
  primer_search_delete(search);
  /* the errors of a tail ending at a position */
  test_eq((int)tail_errors("GGGGGG", 6, "ACGGGGGGT", 8), 0, int, "%d");
  test_eq((int)tail_errors("GGGGGG", 6, "ACGGGGGGT", 9), 1, int, "%d");
  test_eq((int)tail_errors("GGGGGG", 6, "ACGGAGGGT", 8), 1, int, "%d");
  test_eq((int)tail_errors("GGGGGG", 6, "GGG", 3), 3, int, "%d");
  scheme_choice_free(&c);
  delete_logger(logger);
}
