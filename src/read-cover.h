/* (c) 2026 Giorgio Gonnella */

/*
   Where a read lies on its reference, and its category.

   The span of a read on a reference is [start, start + template length)
   on the plus strand, anchored at whichever end mapped (see amplicon.h:
   the mapped probe positions and the kept length), or else where the
   read was placed by a probe from inside it (placement.h), or, for a
   complete amplicon whose probes did not map, the whole segment of its
   length (if its length tells which segment it is: within 80 bases of
   that of a reference, and not as near to a reference of another
   segment). This is an overview from the classification, not a
   base-by-base mapping: it ignores indels.

   The class of a read (COVER_*) is the origin of its 5' (left) end and
   its strand. A read with both ends mapped much farther apart than its
   template length has an internal deletion (COVER_DELETION): it covers
   only its two flanks, up to its junction (junction.h), or, if the
   junction was not located, half of the template at each mapped end.

   The classes are those of the coverage plot of the report
   (counter_coverage()). The category of a read (READ_CATEGORY_*) merges
   them into the four kinds of reads told apart in the output: the tag
   ct of the records, and the files of the --by-category option (see
   trimmer.h).
*/

#ifndef READ_COVER_H
#define READ_COVER_H

#include <stdbool.h>
#include <stddef.h>
#include "amplicon.h"
#include "fasta.h"
#include "reference.h"

/* the classes of a read for the coverage plot: where its amplicon
   started, and on which strand */
enum { COVER_TERMINAL_PLUS,   /* started at a segment end, plus strand */
       COVER_TERMINAL_MINUS,  /* started at a segment end, minus strand */
       COVER_INTERNAL_PLUS,   /* internally primed, plus strand */
       COVER_INTERNAL_MINUS,  /* internally primed, minus strand */
       COVER_NO_PRIMER,       /* no primer at the 5' end (placed by the
                                 other end and the length, or by a probe
                                 from inside the read) */
       COVER_UNMAPPED,        /* a 5' primer, but its probe not located */
       COVER_DELETION,        /* both ends mapped, much farther apart than
                                 the template length: an internal deletion
                                 (e.g. a defective interfering genome) */
       COVER_N_CATEGORIES };

/* the categories of a read in the output */
enum { READ_CATEGORY_FULL,      /* a full template (terminal, any strand) */
       READ_CATEGORY_INTERNAL,  /* internal priming (any strand) */
       READ_CATEGORY_DELETION,  /* an internal deletion */
       READ_CATEGORY_OTHER,     /* no 5' primer, or not located */
       READ_N_CATEGORIES };

/* where a read lies on its reference */
typedef struct {
  bool located;       /* false: the other fields are not set */
  size_t cover;       /* COVER_* */
  size_t ref;         /* the index of the reference */
  /* the span [from, to), 0-based on the plus strand; it is not clamped
     to the reference */
  long from, to;
  /* COVER_DELETION: the flanks are [from, left_to) and [right_from, to) */
  long left_to, right_from;
} read_cover_t;

/*
   Locate the read analyzed as <a> on the references <refs> (NULL: none,
   no read is then located); <set>: the reference set of the analysis,
   whose restrictions to some references and to a segment hold for a read
   placed by its length (NULL: none); false, with <cover>-><located>
   false, if the read cannot be located
*/
bool read_cover_locate(const fasta_t *refs, const reference_set_t *set,
    const amplicon_t *a, read_cover_t *cover);

/* The category (READ_CATEGORY_*) of a read */
size_t read_cover_category(const read_cover_t *cover);

/* The name of a category: full, internal, deletion, other */
const char* read_category_name(size_t category);

/* Unit tests */
void read_cover_test(void);

#endif
