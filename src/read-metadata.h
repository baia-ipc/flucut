/* (c) 2026 Giorgio Gonnella */

/*
   What the titles of the reads say about the sequencing run.

   MinKNOW and Dorado write the metadata of a read in its FASTQ title as
   fields key=value (runid=... flow_cell_id=...
   basecall_model_version_id=...): the flow cell, the run, the sample
   and experiment names given in MinKNOW, the basecalling model (which
   names the pore, the chemistry, the speed and the accuracy, e.g.
   dna_r10.4.1_e8.2_400bps_sup@v4.3.0) and the start time of the read.
   The title of an Illumina read starts with
   <instrument>:<run>:<flow cell>:<lane>:... (Casava 1.8).

   The values are taken from some of the reads (see counter.h); a value
   which differs among them is reported as several (e.g. the reads of
   several runs in one input). Of the start times, the earliest is kept.
*/

#ifndef READ_METADATA_H
#define READ_METADATA_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

/* the values kept */
enum { READ_METADATA_FLOW_CELL, READ_METADATA_RUN, READ_METADATA_MODEL,
       READ_METADATA_SAMPLE, READ_METADATA_EXPERIMENT,
       READ_METADATA_INSTRUMENT, READ_METADATA_START, READ_METADATA_N };

#define READ_METADATA_VALUE_SIZE 80

typedef struct {
  /* ("": not known) */
  char values[READ_METADATA_N][READ_METADATA_VALUE_SIZE];
  bool several[READ_METADATA_N];
} read_metadata_t;

/* Forget all values */
void read_metadata_clear(read_metadata_t *self);

/* Add the values of the title <title> (of length <len>, with or without
   its @) of a read */
void read_metadata_add_title(read_metadata_t *self, const char *title,
    size_t len);

/* Add the values of <other> */
void read_metadata_merge(read_metadata_t *self,
    const read_metadata_t *other);

/* Whether any value is known */
bool read_metadata_known(const read_metadata_t *self);

/* The name of the value <i> for programs (e.g. "flow_cell") and for
   people (e.g. "flow cell") */
const char* read_metadata_key(size_t i);
const char* read_metadata_label(size_t i);

/* The value <i> for people: the value, "several" or NULL (not known) */
const char* read_metadata_text(const read_metadata_t *self, size_t i);

/*
   Save the values to <file>, as lines of text starting with "meta"; and
   load such a line, added to the values: false if <line> is not a valid
   one
*/
void read_metadata_save(const read_metadata_t *self, FILE *file);
bool read_metadata_load_line(read_metadata_t *self, const char *line);

/* Unit tests */
void read_metadata_test(void);

#endif
