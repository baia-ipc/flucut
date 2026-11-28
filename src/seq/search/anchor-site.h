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
   The same for the anchor of <len> bases which begins at the 0-based
   position <pos> of the reference, for a caller which located it
   itself. False, and <self> not found, if the position is outside the
   sequence, if it does not have <len> bases there, or if <len> is 0 or
   more than ANCHOR_MAX_LEN.
*/
bool anchor_at(anchor_t *self, size_t ref, const char *seq, size_t seq_len,
    long pos, size_t len);

/* where an anchor was found in a read */
typedef struct {
  bool found;        /* false: the other fields are not set */
  bool reverse;      /* the read is the reverse complement of the
                        reference, so what precedes the anchor there
                        follows the match here */
  long start, end;   /* 0-based positions of the match in the read, the
                        last one included */
} anchor_hit_t;

/*
   Find the anchor in the <read_len> bases of <read>, in both
   orientations, each with at most <max_error_rate> of its length as
   errors. False, and <hit> not found, if it is not found, or if it is
   found as well in both orientations: a read which matches equally well
   either way does not say which side of the match a window is on.

   One search per read: the windows beside the match are then read with
   anchor_bases_before() and anchor_bases_after(), which do not search
   again.
*/
bool anchor_find(const anchor_t *self, const char *read, size_t read_len,
    double max_error_rate, anchor_hit_t *hit);

/* The bases the anchor matched in the read: its length there, which
   differs from the length of the anchor where the read has an indel in
   it */
static inline size_t anchor_hit_len(const anchor_hit_t *hit) {
  return (size_t)(hit->end - hit->start + 1);
}

/*
   Read the <n> bases of the read which begin <offset> bases from the
   first base of the anchor, counted in the orientation of the reference
   (negative: before it), as the read carries them, into <out> (room for
   n + 1 characters; terminated, upper case). False if the read does not
   reach that far.

   The window is placed by the anchor and read as the read carries it,
   so an insertion or a deletion between the two shifts it: a window
   meant to be read in frame belongs directly beside the anchor.
*/
bool anchor_bases_window(const anchor_hit_t *hit, const char *read,
    size_t read_len, long offset, size_t n, char *out);

/*
   Read the <n> bases which precede the anchor in the orientation of the
   reference, as the read carries them, into <out> (room for n + 1
   characters; terminated, upper case). False if the read does not have
   <n> bases on that side.
*/
bool anchor_bases_before(const anchor_hit_t *hit, const char *read,
    size_t read_len, size_t n, char *out);

/*
   The <n> first bases of the anchor itself, as the read carries them:
   what the read has where the anchor matched. The anchor is conserved,
   so these are a check that it landed where it should, and they show
   the variant of it this sample carries.
*/
bool anchor_bases_at(const anchor_hit_t *hit, const char *read,
    size_t read_len, size_t n, char *out);

/*
   The same for the <n> bases which follow the anchor, <skip> bases
   after its end: with <skip> 0 the ones directly after it, with <skip>
   n the ones after those. Those bases are ordinary sequence beside the
   anchor, which makes them a control for a window read before it.
*/
bool anchor_bases_after(const anchor_hit_t *hit, const char *read,
    size_t read_len, size_t skip, size_t n, char *out);

/* Unit tests */
void anchor_site_test(void);

#endif
