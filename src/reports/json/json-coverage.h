/* (c) 2026 Giorgio Gonnella */

/*
   The coverage by read origin in the JSON summary: for each reference
   with located reads, the depth of each class of reads (read-cover.h)
   at up to JSON_COVERAGE_POINTS positions, evenly spaced along the
   reference. It is the coverage estimated from where the reads are
   located (counter.h), as drawn in the HTML report, so that a program
   can use it without reading the plots.
*/

#ifndef JSON_COVERAGE_H
#define JSON_COVERAGE_H

#include <stdio.h>
#include "counter.h"
#include "fasta.h"

/* the positions at which the depth of a reference is given, at most */
#define JSON_COVERAGE_POINTS 200

/* Whether <c> has a depth on the reference <ref> of <refs> to write */
bool json_coverage_has(const counter_t *c, size_t ref,
    const fasta_t *refs);

/* Write the coverage of the reference <ref> of <refs>, as an object on
   one line: {"reference": ..., "length": ..., "positions": [...],
   "depth": {"terminal_plus": [...], ...}}; the positions are 0-based; a
   class without reads on the reference is left out */
void json_coverage_write(FILE *f, const fasta_t *refs, const counter_t *c,
    size_t ref);

#endif
