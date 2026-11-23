/* (c) 2026 Giorgio Gonnella */

/*
   The reads themselves in the JSON summary: the values of the members
   read_lengths and read_quality (the lengths and the qualities of the
   reads of the input and of the reads written, read-stats.h) and
   sequencing (what the titles of the reads say about the run,
   read-metadata.h), each written on one line.
*/

#ifndef JSON_READS_H
#define JSON_READS_H

#include <stdio.h>
#include "counter.h"

/* {"input": {"reads": ..., "bases": ..., "min": ..., "max": ...,
   "mean": ..., "median": ..., "n50": ...}, "written": {...}} */
void json_reads_lengths(FILE *f, const counter_t *c);

/* {"input": {"reads": ..., "median": ..., "q10": ..., "q20": ...},
   "written": {...}}: the reads with qualities, the median of their
   qualities and the shares of them from Q10 and from Q20; null for a set
   without qualities */
void json_reads_quality(FILE *f, const counter_t *c);

/* {"flow_cell": ..., "run": ..., ...}: each a string, "several" if the
   reads differ, or null if not known */
void json_reads_sequencing(FILE *f, const counter_t *c);

#endif
