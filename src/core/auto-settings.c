/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include "alloc.h"
#include "amplicon.h"
#include "auto-settings.h"
#include "fastq-reader.h"
#include "numutils.h"
#include "read-agreement.h"
#include "reads-reader.h"
#include "reference.h"
#include "segments.h"
#include "strutils.h"
#include "test.h"

void auto_settings_init(auto_settings_t *self, const options_t *opts) {
  assert(self != NULL);
  assert(opts != NULL);
  bool is_auto = opts->max_errors == OPTIONS_MAX_ERRORS_AUTO;
  /* Illumina reads: AUTO_MAX_ERRORS_ILLUMINA, not measured (see
     auto-settings.h) */
  bool automatic = is_auto && !opts->illumina;
  self->max_errors = (report_setting_t){.automatic = automatic,
      .pending = automatic,
      .value = !is_auto ? opts->max_errors :
               opts->illumina ? AUTO_MAX_ERRORS_ILLUMINA : AUTO_MAX_ERRORS_LOW,
      .measured = 0, .n_reads = 0, .file = NULL};
  self->file = NULL;
}

bool auto_settings_pending(const auto_settings_t *self) {
  assert(self != NULL);
  return self->max_errors.pending;
}

double auto_settings_max_errors(double error_rate) {
  return error_rate >= AUTO_MAX_ERRORS_HIGH_RATE ? AUTO_MAX_ERRORS_HIGH :
         AUTO_MAX_ERRORS_LOW;
}

/* add the errors and the aligned read bases of an end found */
/*
   Add the primer found at an end to the errors and the bases of the
   type of influenza that primer belongs to (<primers>: the primers of
   that end of the scheme)
*/
static void add_end(const amplicon_end_t *end, const primer_t *primers,
    size_t *errors, size_t *bases) {
  if (!end->found)
    return;
  virus_type_t virus = primers[end->variant].virus;
  errors[virus] += (size_t)end->hit.n_errors;
  bases[virus] += (size_t)(end->hit.end - end->hit.start + 1);
}

double auto_settings_error_rate(const trimmer_config_t *config,
    const char *file, size_t max_reads, size_t *n_reads, Logger *logger) {
  assert(config != NULL);
  assert(file != NULL);
  amplicon_analyzer_t *a = amplicon_analyzer_new(config->scheme,
      AUTO_MAX_ERRORS_LOW, config->policy, config->window);
  reads_reader_t *reader = reads_reader_new(file, logger);
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  /* the errors and the aligned bases of each type of influenza: a
     scheme can hold the primers of both, and the type of the run is
     not known yet (the references say it only once a sample is read,
     see auto_settings_choose_references). The primers of the other
     type are not in the reads, so the few ends where one of them
     matches best match it with many errors; counting them with the
     rest would raise the rate measured, and the rate decides the
     errors allowed in a primer. The rate returned is therefore that of
     the type whose primers cover most of the read ends. */
  size_t n = 0, got, errors[VIRUS_N_TYPES] = {0};
  size_t bases[VIRUS_N_TYPES] = {0};
  while (n < max_reads &&
         (got = reads_reader_read_batch(reader, &batch, max_reads - n)) > 0) {
    for (size_t i = 0; i < got; i++) {
      const fastq_record_t *rec = batch.records + i;
      const amplicon_t *r = amplicon_analyzer_run(a,
          batch.data.s + rec->seq_from, rec->seq_len);
      add_end(&r->left, config->scheme->left, errors, bases);
      add_end(&r->right, config->scheme->right, errors, bases);
    }
    n += got;
  }
  fastq_batch_free(&batch);
  reads_reader_delete(reader);
  amplicon_analyzer_delete(a);
  if (n_reads != NULL)
    *n_reads = n;
  virus_type_t best = VIRUS_A;
  for (int virus = 1; virus < VIRUS_N_TYPES; virus++)
    if (bases[virus] > bases[best])
      best = (virus_type_t)virus;
  return bases[best] == 0 ? 0.0 :
         (double)errors[best] / (double)bases[best];
}

/* the reference with most votes among the <n> references of the segment
   <s> (of references with as many votes, the first), SIZE_MAX if the
   segment has none; their number and their votes are stored in
   <n_refs> and <total> */
static size_t segment_best(const int *segment, const size_t *votes, size_t n,
    int s, size_t *n_refs, size_t *total) {
  size_t best = SIZE_MAX;
  *n_refs = 0;
  *total = 0;
  for (size_t i = 0; i < n; i++) {
    if (segment[i] != s)
      continue;
    if (best == SIZE_MAX || votes[i] > votes[best])
      best = i;
    (*n_refs)++;
    *total += votes[i];
  }
  return best;
}

/* whether the segment whose leading subtype (or reference, where the
   references have none) has <best_votes> of the <total> votes of the
   segment is clear */
static bool segment_is_clear(size_t best_votes, size_t total) {
  return best_votes >= AUTO_REFERENCE_MIN_VOTES &&
         auto_settings_clear_lead(best_votes, total);
}

/* The subtype with most votes among the references of the segment <s>,
   and its votes (in <lead_votes>): the votes of the references of a
   subtype are counted together, as several references of one subtype
   (its lineages of different decades, hosts and places) share them,
   and one reference of another subtype would otherwise lead with
   fewer. -1 if no reference of the segment has a subtype; the votes of
   the references without one are left out of <lead_votes>. */
static int segment_subtype_lead(const int *segment, const int *subtype,
    const size_t *votes, size_t n, int s, size_t *lead_votes) {
  int lead = -1;
  *lead_votes = 0;
  if (subtype == NULL)
    return -1;
  for (size_t i = 0; i < n; i++) {
    if (segment[i] != s || subtype[i] < 0)
      continue;
    size_t sum = 0;
    for (size_t j = 0; j < n; j++)
      if (segment[j] == s && subtype[j] == subtype[i])
        sum += votes[j];
    if (lead < 0 || sum > *lead_votes) {
      lead = subtype[i];
      *lead_votes = sum;
    }
  }
  return lead;
}

/* the reference with most votes among those of the segment <s> with the
   subtype <sub> (of references with as many votes, the first),
   SIZE_MAX if there is none */
static size_t subtype_best(const int *segment, const int *subtype,
    const size_t *votes, size_t n, int s, int sub) {
  size_t best = SIZE_MAX;
  for (size_t i = 0; i < n; i++) {
    if (segment[i] != s || subtype[i] != sub)
      continue;
    if (best == SIZE_MAX || votes[i] > votes[best])
      best = i;
  }
  return best;
}

/* The reference to use for the segment <s>: the one with most votes
   among those of the subtype with most votes in the segment, or the
   one with most votes where the references have no subtype. The
   subtype decides first, because the lineages of one subtype share the
   votes of their subtype and one reference of another subtype would
   otherwise lead with fewer.

   <lead> and <lead_votes> take the subtype and its votes (the
   reference's own votes where there is no subtype), <total> the votes
   of the segment and <clear> whether the lead is clear (see
   segment_is_clear()). SIZE_MAX if the segment has fewer than two
   references. */
static size_t segment_choice(const int *segment, const int *subtype,
    const size_t *votes, size_t n, int s, int *lead, size_t *lead_votes,
    size_t *total, bool *clear) {
  size_t n_refs;
  size_t best = segment_best(segment, votes, n, s, &n_refs, total);
  *lead = -1;
  *lead_votes = 0;
  *clear = false;
  if (best == SIZE_MAX || n_refs < 2)
    return SIZE_MAX;
  *lead = segment_subtype_lead(segment, subtype, votes, n, s, lead_votes);
  if (*lead < 0) {
    *lead_votes = votes[best];
    *clear = segment_is_clear(votes[best], *total);
    return best;
  }
  *clear = segment_is_clear(*lead_votes, *total);
  return subtype_best(segment, subtype, votes, n, s, *lead);
}

int auto_settings_common_subtype(const int *segment, const int *subtype,
    const size_t *votes, size_t n) {
  if (subtype == NULL)
    return -1;
  /* the subtype of each clear segment and its votes there */
  int chosen[SEGMENTS_N + 1];
  size_t chosen_votes[SEGMENTS_N + 1];
  for (int s = 1; s <= SEGMENTS_N; s++) {
    int lead = -1;
    size_t lead_votes = 0, total;
    bool clear = false;
    size_t pick = segment_choice(segment, subtype, votes, n, s, &lead,
                                 &lead_votes, &total, &clear);
    chosen[s] = pick != SIZE_MAX && clear ? lead : -1;
    chosen_votes[s] = pick != SIZE_MAX && clear ? lead_votes : 0;
  }
  /* the subtype of most clear segments; then the one with most votes in
     them; then the one of the first segment */
  int common = -1;
  size_t common_segments = 0, common_votes = 0;
  for (int s = 1; s <= SEGMENTS_N; s++) {
    if (chosen[s] < 0)
      continue;
    size_t n_segments = 0, n_votes = 0;
    for (int s2 = 1; s2 <= SEGMENTS_N; s2++) {
      if (chosen[s2] == chosen[s]) {
        n_segments++;
        n_votes += chosen_votes[s2];
      }
    }
    if (n_segments > common_segments ||
        (n_segments == common_segments && n_votes > common_votes)) {
      common = chosen[s];
      common_segments = n_segments;
      common_votes = n_votes;
    }
  }
  return common;
}

/* the reference of the subtype <common> to use for the segment <s>,
   whose reference with most votes is <best>: the one with most votes
   among those which are possible choices (all with too few votes, the
   close contenders otherwise); SIZE_MAX if there is none */
static size_t follow_common(const int *segment, const int *subtype,
    const size_t *votes, size_t n, int s, size_t best, int common) {
  /* (the HA and the NA make the subtype: they are chosen by their own
     reads only, never after the other segments) */
  if (common < 0 || s == AUTO_REFERENCE_HA || s == AUTO_REFERENCE_NA)
    return SIZE_MAX;
  bool few = votes[best] < AUTO_REFERENCE_MIN_VOTES;
  size_t chosen = SIZE_MAX;
  for (size_t i = 0; i < n; i++) {
    if (segment[i] != s || subtype[i] != common)
      continue;
    if (!few && (double)votes[i] <
                AUTO_REFERENCE_CLOSE_SHARE * (double)votes[best])
      continue;
    if (chosen == SIZE_MAX || votes[i] > votes[chosen])
      chosen = i;
  }
  return chosen;
}

size_t auto_settings_decide_references(const int *segment,
    const int *subtype, const size_t *votes, size_t n, bool *allowed) {
  size_t undecided = 0;
  int common = auto_settings_common_subtype(segment, subtype, votes, n);
  for (size_t i = 0; i < n; i++)
    allowed[i] = true;
  for (int s = 1; s <= SEGMENTS_N; s++) {
    size_t n_refs, total;
    size_t best = segment_best(segment, votes, n, s, &n_refs, &total);
    if (n_refs < 2)
      continue;
    int lead = -1;
    size_t lead_votes = 0, unused;
    bool clear = false;
    size_t chosen = segment_choice(segment, subtype, votes, n, s, &lead,
                                   &lead_votes, &unused, &clear);
    if (!clear) {
      /* as the clear segments, if possible; otherwise the reference of
         the leading subtype, as more reads would not change the
         shares; with too few votes there is no choice */
      size_t follow = follow_common(segment, subtype, votes, n, s, best,
                                    common);
      if (follow != SIZE_MAX) {
        chosen = follow;
      } else if (lead_votes < AUTO_REFERENCE_MIN_VOTES) {
        undecided++;
        continue;
      }
    }
    for (size_t i = 0; i < n; i++)
      if (segment[i] == s)
        allowed[i] = i == chosen;
  }
  return undecided;
}

bool auto_settings_clear_lead(size_t best_votes, size_t segment_votes) {
  return (double)best_votes >=
         AUTO_REFERENCE_MIN_SHARE * (double)segment_votes;
}

/* add to <votes> those of the probes of the read <seq> (length <len>),
   see auto_settings_choose_references(); the probes which do not vote,
   as several references match them, are counted in <ambiguous>, by the
   segment <segment> of the reference found; the votes of the read are
   added to <agreement> and, by the HA or NA subtype <number> of the
   references, to <subtype_reads> */
static void vote(reference_set_t *refs, const int *segment,
    const int *number, const char *seq, size_t len, size_t *votes,
    size_t *ambiguous, read_agreement_t *agreement,
    subtype_reads_t *subtype_reads) {
  size_t probe_len = reference_set_probe_len(refs);
  if (len < probe_len)
    return;
  size_t room = len - probe_len, last = len;
  size_t read_votes[AUTO_REFERENCE_PROBES], n_read_votes = 0;
  for (size_t j = 0; j < AUTO_REFERENCE_PROBES; j++) {
    size_t offset = room * (2 * j + 1) / (2 * AUTO_REFERENCE_PROBES);
    if (offset == last)
      continue;   /* (a short read: the same probe again) */
    last = offset;
    reference_hit_t hit;
    reference_set_locate_seeded(refs, seq + offset, probe_len, &hit);
    if (!hit.mapped)
      continue;
    if (hit.ambiguous)
      ambiguous[segment[hit.ref_index]]++;
    else {
      votes[hit.ref_index]++;
      read_votes[n_read_votes++] = hit.ref_index;
    }
  }
  read_agreement_add(agreement, segment, read_votes, n_read_votes);
  subtype_reads_add(subtype_reads, segment, number, read_votes,
                    n_read_votes);
}

virus_set_t auto_settings_reference_types(const trimmer_config_t *config,
    const bool *allowed) {
  assert(config != NULL);
  if (config->index == NULL)
    return VIRUS_SET_BOTH;
  const fasta_t *fa = reference_index_fasta(config->index);
  virus_set_t types = 0;
  for (size_t i = 0; i < fa->n_records; i++) {
    if (allowed != NULL && !allowed[i])
      continue;
    const char *description = fa->records[i].description;
    types |= VIRUS_SET_OF(description == NULL ? VIRUS_A :
                          virus_type_of_description(description));
  }
  return types == 0 ? VIRUS_SET_BOTH : types;
}

/* the description of a reference of the segment <s>, for the numbering
   of its type in the messages (NULL if there is none) */
static const char* segment_description(const fasta_t *fa,
    const int *segment, int s) {
  for (size_t i = 0; i < fa->n_records; i++)
    if (segment[i] == s)
      return fa->records[i].description;
  return NULL;
}

/* report the references chosen for <file> */
static void log_references(const trimmer_config_t *config, const char *file,
    const int *segment, const int *subtype, const size_t *votes,
    const size_t *ambiguous, const bool *allowed, size_t n_reads,
    Logger *logger) {
  const fasta_t *fa = reference_index_fasta(config->index);
  int common = auto_settings_common_subtype(segment, subtype, votes,
                                            fa->n_records);
  for (int s = 1; s <= SEGMENTS_N; s++) {
    size_t n_refs, total, n_allowed = 0, chosen = 0;
    size_t best = segment_best(segment, votes, fa->n_records, s, &n_refs,
                               &total);
    if (n_refs < 2)
      continue;
    for (size_t i = 0; i < fa->n_records; i++) {
      if (segment[i] == s && allowed[i]) {
        chosen = i;
        n_allowed++;
      }
    }
    char name[16];
    segments_name_in(s, segment_description(fa, segment, s), name,
                     sizeof(name));
    for (size_t i = 0; i < fa->n_records; i++)
      if (segment[i] == s)
        lprintf(logger, LOG_DEBUG, "votes of %s for %s: %s %zu\n", name,
                file, fa->records[i].name, votes[i]);
    if (n_allowed != 1) {
      lprintf(logger, LOG_INFO, "reference of %s for %s: not chosen, too "
              "few probes (%zu) in the first %zu reads (%zu more match "
              "several references): all %zu are used\n", name, file, total,
              n_reads, ambiguous[s], n_refs);
      continue;
    }
    /* how it was chosen: by the votes of its subtype in this segment,
       or as the clear segments */
    int lead = -1;
    size_t lead_votes = 0, of_segment;
    bool clear = false;
    segment_choice(segment, subtype, votes, fa->n_records, s, &lead,
                   &lead_votes, &of_segment, &clear);
    const char *how = "";
    char subtype_votes[128] = "";
    if (lead >= 0 && lead_votes > votes[chosen])
      snprintf(subtype_votes, sizeof(subtype_votes),
               ", %zu for its subtype", lead_votes);
    if (!clear) {
      if (common >= 0 && subtype[chosen] == common)
        how = votes[best] < AUTO_REFERENCE_MIN_VOTES ?
              "; too few probes: the subtype of the clear segments" :
              "; no clear lead: the subtype of the clear segments";
      else
        how = "; about as close to other references";
    }
    lprintf(logger, LOG_INFO, "reference of %s for %s: %s (%zu of %zu "
            "probes of the first %zu reads%s; %zu more match several "
            "references)%s\n", name, file, fa->records[chosen].name,
            votes[chosen], total, n_reads, subtype_votes, ambiguous[s],
            how);
  }
}

/*
   The type of the sample from the votes: the type of the references in
   which the probes were located, when it has at least
   AUTO_VIRUS_SHARE of the votes; both types otherwise (a sample of
   both, or too few votes to tell), so that the primers of both are
   searched when it is not clear which one it is.
*/
static virus_set_t votes_virus(const fasta_t *fa, const size_t *votes,
    size_t n) {
  size_t by_type[VIRUS_N_TYPES] = {0}, total = 0;
  for (size_t i = 0; i < n; i++) {
    if (votes[i] == 0)
      continue;
    const char *description = fa->records[i].description;
    virus_type_t type = description == NULL ? VIRUS_A :
                        virus_type_of_description(description);
    by_type[type] += votes[i];
    total += votes[i];
  }
  if (total == 0)
    return VIRUS_SET_BOTH;
  for (int type = 0; type < VIRUS_N_TYPES; type++)
    if ((double)by_type[type] >= AUTO_VIRUS_SHARE * (double)total)
      return VIRUS_SET_OF((virus_type_t)type);
  return VIRUS_SET_BOTH;
}

bool* auto_settings_choose_references(const trimmer_config_t *config,
    const char *file, virus_set_t *types, subtype_reads_t *subtype_reads,
    Logger *logger) {
  assert(config != NULL);
  assert(file != NULL);
  /* the types of all the references, until the votes say which one the
     sample is */
  if (types != NULL)
    *types = auto_settings_reference_types(config, NULL);
  subtype_reads_t of_subtypes;
  subtype_reads_clear(&of_subtypes);
  if (subtype_reads != NULL)
    *subtype_reads = of_subtypes;
  if (config->index == NULL)
    return NULL;
  const fasta_t *fa = reference_index_fasta(config->index);
  size_t n_refs = fa->n_records;
  /* a segment with several references is needed for a choice */
  size_t per_segment[SEGMENTS_N + 1] = {0};
  bool several = false;
  int *segment, *subtype;
  XMALLOC(segment, sizeof(int) * (n_refs > 0 ? n_refs : 1));
  XMALLOC(subtype, sizeof(int) * (n_refs > 0 ? n_refs : 1));
  for (size_t i = 0; i < n_refs; i++) {
    subtype[i] = reference_index_subtype(config->index, i);
    segment[i] = reference_index_segment(config->index, i);
    if (segment[i] > 0 && ++per_segment[segment[i]] > 1)
      several = true;
  }
  if (!several) {
    free(segment);
    free(subtype);
    return NULL;
  }
  reference_set_t *refs = reference_set_new(config->index, config->probe_len,
      config->probe_error_rate);
  size_t *votes, ambiguous[SEGMENTS_N + 1] = {0};
  bool *allowed;
  XCALLOC(votes, n_refs, sizeof(size_t));
  XMALLOC(allowed, sizeof(bool) * n_refs);
  read_agreement_t agreement;
  read_agreement_init(&agreement, n_refs);
  int *number = subtype_reads_numbers(fa, segment);
  of_subtypes.samples = 1;
  reads_reader_t *reader = reads_reader_new(file, logger);
  fastq_batch_t batch;
  fastq_batch_init(&batch);
  size_t n = 0, undecided = 0;
  const size_t stages[] = {AUTO_REFERENCE_READS, AUTO_REFERENCE_MORE_READS};
  for (size_t stage = 0; stage < 2; stage++) {
    size_t got = 0;
    while (n < stages[stage] &&
           (got = reads_reader_read_batch(reader, &batch,
                                          stages[stage] - n)) > 0) {
      for (size_t i = 0; i < got; i++) {
        const fastq_record_t *rec = batch.records + i;
        vote(refs, segment, number, batch.data.s + rec->seq_from,
             rec->seq_len, votes, ambiguous, &agreement, &of_subtypes);
      }
      n += got;
    }
    undecided = auto_settings_decide_references(segment, subtype, votes,
                                                n_refs, allowed);
    /* decided, or no more reads */
    if (undecided == 0 || n < stages[stage])
      break;
  }
  fastq_batch_free(&batch);
  reads_reader_delete(reader);
  reference_set_delete(refs);
  if (types != NULL)
    *types = votes_virus(fa, votes, n_refs);
  log_references(config, file, segment, subtype, votes, ambiguous, allowed, n,
                 logger);
  read_agreement_log(&agreement, fa, segment, votes, file, logger);
  read_agreement_free(&agreement);
  subtype_reads_log(&of_subtypes, file, logger);
  if (subtype_reads != NULL)
    *subtype_reads = of_subtypes;
  free(number);
  bool all = true;
  for (size_t i = 0; i < n_refs; i++)
    if (!allowed[i])
      all = false;
  free(votes);
  free(segment);
  free(subtype);
  if (all) {
    free(allowed);
    return NULL;
  }
  return allowed;
}

void auto_settings_decide(auto_settings_t *self, trimmer_config_t *config,
    const char *file, Logger *logger) {
  assert(self != NULL);
  assert(config != NULL);
  assert(file != NULL);
  if (!self->max_errors.pending)
    return;
  report_setting_t *s = &self->max_errors;
  s->measured = auto_settings_error_rate(config, file,
      AUTO_MAX_ERRORS_READS, &s->n_reads, logger);
  s->value = auto_settings_max_errors(s->measured);
  config->max_errors = s->value;
  free(self->file);
  self->file = xstrdup(file);
  s->file = self->file;
  s->pending = false;
  lprintf(logger, LOG_INFO, "Error rate of the primer alignments in the "
          "first %zu reads of %s: %s; errors in a primer: %g%% of its "
          "length (-e auto)\n", s->n_reads, file,
          fraction_percent_text(s->measured).s, 100 * s->value);
}

void auto_settings_free(auto_settings_t *self) {
  if (self == NULL)
    return;
  free(self->file);
  self->file = NULL;
  self->max_errors.file = NULL;
}

void auto_settings_test(void) {
  /* the choice of the reference of each segment: segment 4 has three
     references, segment 6 two, segment 1 one; a reference without a
     segment */
  {
    const int segment[] = {4, 4, 4, 6, 6, 1, 0};
    bool allowed[7];
    /* a clear winner in each segment */
    const size_t clear[] = {3, 90, 7, 40, 0, 0, 0};
    test_eq_size(auto_settings_decide_references(segment, NULL, clear,
                                                 7, allowed), 0);
    test_assert(!allowed[0] && allowed[1] && !allowed[2]);
    test_assert(allowed[3] && !allowed[4]);
    test_assert(allowed[5] && allowed[6]);
    /* too few votes in segment 6: it keeps its references */
    const size_t few[] = {3, 90, 7, AUTO_REFERENCE_MIN_VOTES - 1, 0, 0, 0};
    test_eq_size(auto_settings_decide_references(segment, NULL, few, 7,
                                                 allowed),
                 1);
    test_assert(!allowed[0] && allowed[1] && !allowed[2]);
    test_assert(allowed[3] && allowed[4]);
    /* no clear lead in segment 4 (a tie: the first is chosen), no vote
       at all in segment 6 */
    const size_t tie[] = {50, 50, 0, 0, 0, 0, 0};
    test_eq_size(auto_settings_decide_references(segment, NULL, tie, 7,
                                                 allowed),
                 1);
    test_assert(allowed[0] && !allowed[1] && !allowed[2]);
    test_assert(allowed[3] && allowed[4]);
    test_assert(!auto_settings_clear_lead(50, 100));
    test_assert(auto_settings_clear_lead(60, 100));
    /* just the share needed */
    const size_t share[] = {60, 40, 0, 30, 20, 0, 0};
    test_eq_size(auto_settings_decide_references(segment, NULL, share,
                                                 7, allowed), 0);
    test_assert(allowed[0] && !allowed[1] && !allowed[2]);
    test_assert(allowed[3] && !allowed[4]);
  }
  /* the votes of the references of a subtype are counted together:
     segment 4 has four references of the subtype 0 and one of the
     subtype 1, which leads among the references but not among the
     subtypes */
  {
    const int segment[] = {4, 4, 4, 4, 4};
    const int subtype[] = {0, 0, 0, 0, 1};
    const size_t votes[] = {95, 61, 15, 7, 153};
    bool allowed[5];
    test_eq_size(auto_settings_decide_references(segment, subtype, votes,
                                                 5, allowed), 0);
    /* the subtype 0 has 178 votes of the 331 of the segment: it is
       clear, and its reference with most votes is chosen */
    test_assert(allowed[0] && !allowed[1] && !allowed[2] && !allowed[3]);
    test_assert(!allowed[4]);
    /* with the votes of the subtype 1 above the share, it is chosen */
    const size_t other[] = {30, 20, 10, 5, 200};
    test_eq_size(auto_settings_decide_references(segment, subtype, other,
                                                 5, allowed), 0);
    test_assert(!allowed[0] && allowed[4]);
    /* no subtype above the share: the reference with most votes */
    const size_t split[] = {60, 40, 0, 0, 110};
    test_eq_size(auto_settings_decide_references(segment, subtype, split,
                                                 5, allowed), 0);
    test_assert(!allowed[0] && allowed[4]);
  }
  /* segments without a clear choice follow the clear ones: three
     subtypes (0, 1, 2) with a reference each in the segments 4, 6, 7
     and 8 */
  {
    const int segment[] = {4, 4, 4, 6, 6, 6, 7, 7, 7, 8, 8, 8};
    const int subtype[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
    bool allowed[12];
    /* segments 4 and 6 are clear for subtype 1; segment 7 has too few
       votes; in segment 8 subtype 2 leads, subtype 1 is close */
    const size_t follow[] = {5, 90, 5, 0, 40, 2, 3, 0, 1, 0, 40, 50};
    test_eq(auto_settings_common_subtype(segment, subtype, follow, 12), 1,
            int, "%d");
    test_eq_size(auto_settings_decide_references(segment, subtype, follow,
                                                 12, allowed), 0);
    test_assert(!allowed[0] && allowed[1] && !allowed[2]);
    test_assert(!allowed[3] && allowed[4] && !allowed[5]);
    test_assert(!allowed[6] && allowed[7] && !allowed[8]);
    test_assert(!allowed[9] && allowed[10] && !allowed[11]);
    /* in segment 8 the reference of the common subtype is not close:
       the one with most votes is chosen */
    const size_t far[] = {5, 90, 5, 0, 40, 2, 3, 0, 1, 0, 24, 50};
    test_eq_size(auto_settings_decide_references(segment, subtype, far, 12,
                                                 allowed), 0);
    test_assert(!allowed[9] && !allowed[10] && allowed[11]);
    /* the clear segments disagree: the subtype of most of them, then
       the one with most votes */
    const size_t most[] = {90, 5, 5, 0, 40, 2, 0, 30, 0, 0, 0, 0};
    test_eq(auto_settings_common_subtype(segment, subtype, most, 12), 1,
            int, "%d");
    const size_t votes_tie[] = {90, 5, 5, 0, 40, 2, 0, 0, 0, 0, 0, 0};
    test_eq(auto_settings_common_subtype(segment, subtype, votes_tie, 12),
            0, int, "%d");
    /* no clear segment: as without subtypes */
    const size_t none[] = {50, 50, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0};
    test_eq(auto_settings_common_subtype(segment, subtype, none, 12), -1,
            int, "%d");
    test_eq_size(auto_settings_decide_references(segment, subtype, none, 12,
                                                 allowed), 3);
    test_assert(allowed[0] && !allowed[1] && !allowed[2]);
    test_assert(allowed[3] && allowed[4] && allowed[5]);
    /* the HA and the NA never follow: with the other segments clear
       for subtype 1 and too few votes of their own, they keep all their
       references */
    const int with_internal[] = {4, 4, 4, 6, 6, 6, 7, 7, 7, 8, 8, 8};
    const size_t own[] = {1, 0, 4, 0, 0, 0, 3, 90, 1, 0, 40, 2};
    test_eq(auto_settings_common_subtype(with_internal, subtype, own, 12),
            1, int, "%d");
    test_eq_size(auto_settings_decide_references(with_internal, subtype,
                                                 own, 12, allowed), 2);
    test_assert(allowed[0] && allowed[1] && allowed[2]);
    test_assert(allowed[3] && allowed[4] && allowed[5]);
    test_assert(!allowed[6] && allowed[7] && !allowed[8]);
    /* ... and with votes but no clear lead, they take their leader, not
       the close reference of the common subtype */
    const size_t lead[] = {0, 40, 50, 0, 0, 0, 3, 90, 1, 0, 40, 2};
    auto_settings_decide_references(with_internal, subtype, lead, 12,
                                    allowed);
    test_assert(!allowed[0] && !allowed[1] && allowed[2]);
    /* references without a subtype do not follow */
    const int unknown[] = {0, 1, 2, 0, 1, 2, -1, -1, -1, 0, 1, 2};
    test_eq_size(auto_settings_decide_references(segment, unknown, follow,
                                                 12, allowed), 1);
    test_assert(allowed[6] && allowed[7] && allowed[8]);
  }
  test_eq(auto_settings_max_errors(0.01), AUTO_MAX_ERRORS_LOW, double,
          "%g");
  test_eq(auto_settings_max_errors(AUTO_MAX_ERRORS_HIGH_RATE),
          AUTO_MAX_ERRORS_HIGH, double, "%g");
  test_eq(auto_settings_max_errors(0.1), AUTO_MAX_ERRORS_HIGH, double,
          "%g");
  /* -e auto: measured on Nanopore reads, not on Illumina reads */
  {
    options_t opts = {.max_errors = OPTIONS_MAX_ERRORS_AUTO};
    auto_settings_t autos;
    auto_settings_init(&autos, &opts);
    test_assert(autos.max_errors.automatic && auto_settings_pending(&autos));
    auto_settings_free(&autos);
    opts.illumina = true;
    auto_settings_init(&autos, &opts);
    test_assert(!autos.max_errors.automatic &&
                !auto_settings_pending(&autos));
    test_eq(autos.max_errors.value, AUTO_MAX_ERRORS_ILLUMINA, double, "%g");
    auto_settings_free(&autos);
    opts.max_errors = 0.1;
    auto_settings_init(&autos, &opts);
    test_eq(autos.max_errors.value, 0.1, double, "%g");
    auto_settings_free(&autos);
  }
}
