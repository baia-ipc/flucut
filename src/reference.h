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

   This needs references, and does not guess: without them flucut does not
   run. The references are whole segments (or at least enough of each
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

typedef struct reference_set_t reference_set_t;

/*
   Create a reference set from the sequences <fasta> (which must exist as
   long as the set is used).

   <probe_len>: the number of bases after the 5' primer used to locate the
   read in the references; <max_error_rate>: the fraction of the probe
   allowed to be errors (the probe is matched with at most
   probe_len * max_error_rate errors).
*/
reference_set_t* reference_set_new(const fasta_t *fasta, size_t probe_len,
    double max_error_rate);

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
