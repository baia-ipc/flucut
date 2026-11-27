/* (c) 2026 Giorgio Gonnella */

/*
   A short conserved stretch of a reference, used to read what a read
   carries at a position inside a segment without mapping the read.

   flucut locates the ends of the reads, not their interior: the probe
   after a primer is mapped, the rest of the read is not. To report what
   the reads carry at an interior position which matters - the cleavage
   site of the haemagglutinin (cleavage-site.h), and in the same way any
   other position worth watching - that position has to be found inside
   each read itself.

   It is found by a short stretch of the reference beside it, the
   anchor. An anchor is useful when it is

   - conserved, so that it is in the reads of any strain, and can be
     taken from any reference of the segment;
   - outside the stretch to be read, so that an error in that stretch
     does not move it;
   - immediately beside it, so that an insertion or a deletion between
     the two cannot shift what is read out of frame. The bases reported
     are then simply the bases of the read next to the match, as the
     read carries them, with no alignment of the read to the reference
     and no assembly.

   The anchor is searched with the aligner used for the primers (edlib,
   infix, the IUPAC codes equal to their bases), so that a read with
   errors still finds it, and in both orientations, so that the strand
   of the read does not have to be known. A read in which the two
   orientations match equally well is left out: it does not say which
   side of the match the position is on. The read is searched as it is,
   as the primers are, so a read must be in upper case, as the
   basecallers write it; the bases reported are upper-cased.
*/

#ifndef ANCHOR_SITE_H
#define ANCHOR_SITE_H

#include <stdbool.h>
#include <stddef.h>

/* the longest anchor (a few codons are enough for a conserved peptide) */
#define ANCHOR_MAX_LEN 48

typedef struct {
  bool found;        /* false: the other fields are not set */
  size_t ref;        /* the reference the anchor was taken from */
  long ref_pos;      /* 0-based position of its first base there */
  size_t len;
  char seq[ANCHOR_MAX_LEN + 1];
  char rev[ANCHOR_MAX_LEN + 1];   /* its reverse complement */
} anchor_t;

/*
   Build in <self> the anchor of <len> bases which begins at the first
   codon of <peptide> in the reference <ref>, whose sequence is the
   <seq_len> bases of <seq> (see codons_find_peptide()). False, and
   <self> not found, if the peptide is not in the sequence, if the
   sequence does not have <len> bases there, or if <len> is 0 or more
   than ANCHOR_MAX_LEN.
*/
bool anchor_of_peptide(anchor_t *self, size_t ref, const char *seq,
    size_t seq_len, const char *peptide, size_t len);

/*
   Read, out of the <read_len> bases of <read>, the <n> bases which
   precede the anchor in the orientation of the reference, into <out>
   (which must have room for n + 1 characters; the result is terminated
   and in upper case).

   The anchor is searched in both orientations, each with at most
   <max_error_rate> of its length as errors. False, and <out> not
   written, if it is not found, if it is found as well in both
   orientations, or if the read does not have <n> bases before it.
*/
bool anchor_bases_before(const anchor_t *self, const char *read,
    size_t read_len, double max_error_rate, size_t n, char *out);

/* Unit tests */
void anchor_site_test(void);

#endif
