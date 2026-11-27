/* (c) 2026 Giorgio Gonnella */

/*
   Reads which contain several amplicons, and the parts they are split
   into.

   A read is expected to be one amplicon (amplicon.h). Some reads hold
   two or more, joined end to end: molecules ligated to each other
   during the library preparation, or sequenced one after the other
   without the read being split. The primers of the inner amplicons are
   then inside the read, where the trimming of the ends does not reach
   them, and their bases are taken for template: at an internal priming
   site they carry the artefact into the consensus, as the primer of an
   internally primed read does when it is not removed.

   Such a read is split at its junctions, and each part is then analyzed
   and trimmed as a read of its own.

   The junctions
   -------------

   The conserved termini of the primers of the scheme are searched in
   the whole read, on both strands (terminus-matches.h): a match is
   either the start of an amplicon or the end of one.

   A match inside a read is not always a junction: a segment can contain
   a sequence similar to a primer (the internal priming sites), and a
   read of a complete amplicon runs through it. A junction is accepted
   only on one of these grounds:

   - a pair: the end of an amplicon followed, within CHIMERA_MAX_GAP
     bases, by the start of another (between them, at most the tails,
     adapters and barcodes of the two molecules), the two matches having
     at most CHIMERA_MAX_PAIR_ERRORS errors together (two matches next
     to each other by chance are otherwise too common). The read is cut
     after the first match and before the second; the bases between them
     are dropped. With references, a pair is left alone if the read runs
     through it: the bases before it and those after it are located next
     to each other in a reference, as below;
   - with references, a single match more than <window> bases from the
     end of the read it would belong to, when the read does not run
     through it: the bases on the amplicon side of the match are located
     in the references, and those on its other side are not located next
     to them (on the same reference and strand, at the distance they
     have in the read, give or take CHIMERA_CONTINUITY_TOLERANCE, or
     farther apart in the same order: a read with an internal deletion,
     whose junction can resemble a primer, is one template), nor
     are the bases further away (CHIMERA_OUTER_PROBES probes are tried,
     the tolerance growing with the distance): a segment can differ
     from the references next to a site. The read is cut on the outer
     side of the primer of the match: a primer with a tail is searched
     around its terminus, so that the tail goes with its primer. The
     amplicon side must hold a whole probe.

   Without references only the pairs are found. A match which runs on in
   the same segment on both sides is template, and is left alone.
*/

#ifndef CHIMERA_H
#define CHIMERA_H

#include <stddef.h>
#include "reference.h"
#include "scheme.h"
#include "virus.h"

/* the greatest number of parts a read is split into (the last part is
   whatever is left of the read) */
#define CHIMERA_MAX_PARTS 8

/* the greatest number of bases between the two matches of a pair */
#define CHIMERA_MAX_GAP 100

/* the most errors allowed in the two matches of a pair together */
#define CHIMERA_MAX_PAIR_ERRORS 2

/* how much the distance of the two sides of a match on the reference may
   differ from their distance in the read, for the read to be taken as
   running through the match (insertions and deletions) */
#define CHIMERA_CONTINUITY_TOLERANCE 10

/* the probes tried on the outer side of a single match, before the read
   is taken as stopping there */
#define CHIMERA_OUTER_PROBES 3

/* a part of a read: 0-based, inclusive */
typedef struct {
  long from, to;
} chimera_part_t;

typedef struct chimera_finder_t chimera_finder_t;

/*
   Create a finder (it contains buffers: use one for each thread);
   <scheme> must exist as long as the finder is used. <max_errors>: the
   errors allowed in a primer, a number of bases, or, if smaller than
   1, a fraction of its length (see terminus-matches.h). <window>: see
   above (that of the primers of the ends, see amplicon.h).
*/
chimera_finder_t* chimera_finder_new(const scheme_t *scheme,
    double max_errors, long window);

/* Use the references <refs> (see above); <refs> must exist as long as
   the finder is used, and belongs to the thread of the finder */
/*
   Search only the termini of the primers of the types <types> (see
   amplicon_analyzer_restrict_virus)
*/
void chimera_finder_restrict_virus(chimera_finder_t *self,
    virus_set_t types);

void chimera_finder_use_references(chimera_finder_t *self,
    reference_set_t *refs);

/*
   Split the read <seq> (length <len>): the parts are stored in <parts>
   (room for CHIMERA_MAX_PARTS), in the order of the read, and their
   number is returned (1, the whole read, if no junction is found; 0 for
   an empty read).
*/
size_t chimera_finder_run(chimera_finder_t *self, const char *seq,
    size_t len, chimera_part_t *parts);

/* Delete a finder (NULL is accepted) */
void chimera_finder_delete(chimera_finder_t *self);

/* Unit tests */
void chimera_test(void);

#endif
