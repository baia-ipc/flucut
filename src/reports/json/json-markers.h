/* (c) 2026 Giorgio Gonnella */

/*
   The markers in the JSON summary: for each marker the reads gave
   (markers.h), what they carry at it and how many of them carry each
   value, so that a program can see a residue present in a part of the
   reads only, which a consensus would not show.
*/

#ifndef JSON_MARKERS_H
#define JSON_MARKERS_H

#include "counter.h"
#include "json-writer.h"

/*
   Write the member "markers": the markers looked for in the reads of
   <c>, in the order of the table, each as an object; an empty array
   where none was (no reference, or no anchor in the references)
*/
void json_markers_write(json_writer_t *w, const counter_t *c);

#endif
