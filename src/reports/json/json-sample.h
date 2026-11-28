/* (c) 2026 Giorgio Gonnella */

/*
   The sample itself in the JSON summary: the member "sample" of a
   section, with its state, its subtype, what deserves attention in it
   and its segments (qc-summary.h, qc-flags.h); null without references.
*/

#ifndef JSON_SAMPLE_H
#define JSON_SAMPLE_H

#include "counter.h"
#include "json-writer.h"
#include "report.h"

/* Write the member "sample" of the counts <c> to the open object */
void json_sample_write(json_writer_t *w, const report_t *r,
    const counter_t *c);

#endif
