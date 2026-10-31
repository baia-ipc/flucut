/* (c) 2026 Giorgio Gonnella */

/*
   Classification of a read against the reference segments, to tell a read
   which starts at the conserved terminus of a segment (a normal amplicon)
   from one which starts at an interior position of a segment (internal
   priming, see docs/internal-priming.md).

   The idea
   --------

   A read of a normal amplicon begins with the universal primer, and the
   sequence after the primer is the start of the segment: it aligns to the
   reference at the very beginning of the segment (the conserved terminus
   is the first bases of the segment).

   A read started at an internal site also begins with the primer, but the
   sequence after the primer is taken from the middle of the segment: it
   aligns to the reference at an interior position. The primer-derived
   bases at the 5' end of such a read are not the template there, and are
   the ones which carry the artefact into the consensus.

   So the sequence of the read after its 5' primer (a probe) is located in
   the references: if it maps near the start of a segment, the read is a
   normal amplicon; if it maps to an interior position, the read is
   internally primed, and the primer-derived prefix is trimmed (only that
   prefix: the rest of the read is genuine interior template and is kept,
   see trimmer).

   The classification decides the trimming only when the conserved
   termini are kept (TRIM_TERMINUS, see amplicon.h), which then requires
   the references; with the default trimming, which removes every
   primer, it measures internal priming for the report. The references
   are whole segments (or at least enough of each
   segment to contain both the terminal and the interior mapping
   positions); scripts/fetch_references.py downloads them.

   The same procedure runs against the references of all eight segments,
   which also tells, from the data, in which segments internal priming
   occurs (the artefact which motivated flucut is in HA; whether it occurs
   elsewhere is measured, not assumed).
*/

#ifndef REFERENCE_H
#define REFERENCE_H

#include <stdbool.h>
#include <stddef.h>
#include "fasta.h"
#include "logger.h"

/* a probe mapped to the references */
typedef struct {
  bool mapped;
  size_t ref_index;    /* the reference the probe mapped to */
  long ref_len;        /* its length */
  /* the probe mapped to the minus strand of the reference (its reverse
     complement mapped to the plus strand): the read runs from the end of
     the segment towards its start */
  bool reverse;
  /* the position in the reference (on its plus strand) of the first base
     of the probe, i.e. of the first base of the read after its 5' primer;
     for a reverse match, this is the last position of the alignment on
     the plus strand */
  long ref_start;
  long n_errors;
} reference_hit_t;

/* probes shorter than this are not located (too short to be placed
   reliably in a segment); reference_set_locate() then reports no hit */
#define REFERENCE_MIN_PROBE_LEN 16

/* default length of the probe */
#define REFERENCE_DEFAULT_PROBE_LEN 32

/* default fraction of the probe allowed to be errors (Nanopore reads) */
#define REFERENCE_DEFAULT_MAX_ERROR_RATE 0.2

/* default terminal window (see reference_hit_is_terminal()): generous,
   as the references downloaded from the databases may lack some of the
   bases of the termini, while the internal sites to detect are hundreds
   of bases inside the segment */
#define REFERENCE_DEFAULT_TERMINAL_WINDOW 60

/*
   The k-mer index of the references
   ---------------------------------

   Aligning each probe to every reference costs time in proportion to the
   number of references. The references are therefore indexed by their
   k-mers (kmer-index.h of libcbaia, k = REFERENCE_KMER_LEN): the k-mers
   of a probe, and of its reverse complement, give candidate diagonals,
   and only the windows around them are aligned, with the same rules as
   the exhaustive search. k-mers occurring more than
   REFERENCE_MAX_KMER_OCC times (low complexity) are not used as seeds. A
   probe for which the k-mers give no match (e.g. an error in each k-mer,
   or a sequence absent from the references) is searched exhaustively:
   the index thus only speeds up the probes it finds. It can differ from
   the exhaustive search only where a better match than the one found
   has no exact k-mer. The exhaustive search can also be chosen for all
   probes, to check the index; or it can be left out (seeds only), for
   reads with few errors (Illumina), whose probes, if they belong to a
   reference, share exact k-mers with it, while many probes belong to no
   reference at all (e.g. reads of the host), each of which would be
   aligned to every reference.

   The index is built once and shared, read-only, by the searchers of all
   threads.
*/

#define REFERENCE_KMER_LEN 10
#define REFERENCE_MAX_KMER_OCC 4096

typedef struct reference_index_t reference_index_t;

/* Index the sequences <fasta> (which must exist as long as the index is
   used) */
reference_index_t* reference_index_new(const fasta_t *fasta);

/* The sequences of the index */
const fasta_t* reference_index_fasta(const reference_index_t *self);

/* Delete (NULL is accepted) */
void reference_index_delete(reference_index_t *self);

/* The subtype of the reference <ref> (an id, -1 if its description has
   none), its number of distinct subtypes, and the name of a subtype id
   (NULL if out of range); the subtype is parsed from the FASTA
   description (e.g. the "(H5N1)" of an influenza record) */
int reference_index_subtype(const reference_index_t *self, size_t ref);
size_t reference_index_n_subtypes(const reference_index_t *self);
const char* reference_index_subtype_name(const reference_index_t *self,
    int subtype);

/* a searcher of probes in the references */
typedef struct reference_set_t reference_set_t;

/*
   Create a searcher for the indexed references <index> (which must exist
   as long as the searcher is used).

   <probe_len>: the number of bases after the 5' primer used to locate the
   read in the references; <max_error_rate>: the fraction of the probe
   allowed to be errors (the probe is matched with at most
   probe_len * max_error_rate errors).
*/
reference_set_t* reference_set_new(const reference_index_t *index,
    size_t probe_len, double max_error_rate);

/* Align each probe to every reference, without the index (slower; the
   result is the one the index approximates) */
void reference_set_set_exhaustive(reference_set_t *self, bool exhaustive);

/* Consider only the references of the subtype <subtype> (an id, see
   reference_index_subtype); -1 (the default) considers all. Each read
   should be classified against the subtype it belongs to, as the
   references of other subtypes match it spuriously (their untranslated
   regions are conserved, their coding regions diverge) */
void reference_set_restrict(reference_set_t *self, int subtype);

/* Search the probes by their k-mers only, without the exhaustive search
   of the probes they do not locate (see above) */
void reference_set_set_seeds_only(reference_set_t *self, bool seeds_only);

/*
   Locate <probe> (length <probe_len>; the start of the read after its 5'
   primer) in the references, on both strands; only the first
   min(<probe_len>, probe length of the set) bases are used, and a probe
   shorter than REFERENCE_MIN_PROBE_LEN is not located.

   The best match (fewest errors) is stored in <hit>. Among equally good
   matches, the one nearest to the terminus from which the read would
   start (the start of the segment for a forward match, its end for a
   reverse one) is chosen: an ambiguous read is classified as terminal,
   and thus kept intact. The searcher contains buffers: use one per
   thread.
*/
void reference_set_locate(reference_set_t *self, const char *probe,
    size_t probe_len, reference_hit_t *hit);

/*
   Whether the mapped probe <hit> is at a terminus of the segment (a
   normal amplicon) rather than at an interior position (internal
   priming): within <terminal_window> bases of the start of the segment
   for a forward match, of its end for a reverse match. The hit must be
   mapped.
*/
bool reference_hit_is_terminal(const reference_hit_t *hit,
    long terminal_window);

/* The number of references */
size_t reference_set_size(const reference_set_t *self);

/* The probe length of the set */
size_t reference_set_probe_len(const reference_set_t *self);

/* Delete (NULL is accepted) */
void reference_set_delete(reference_set_t *self);

/* Unit tests */
void reference_test(void);

#endif
