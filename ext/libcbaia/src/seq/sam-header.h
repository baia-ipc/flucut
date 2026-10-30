/* (c) 2026 Giorgio Gonnella */

/*
   SAM header text of BAM files written from BAM input (e.g. the unaligned
   BAM of Dorado), keeping the provenance of the reads:
   - the read groups (@RG) of the inputs, so that the records keep the RG
     tag of their input; a program which distributes the reads to several
     destinations (e.g. a demultiplexer) can rename them per destination;
   - the programs (@PG) of the inputs, followed by the @PG line of the
     program, linked to the last one of the inputs (PP), as samtools does.

   The header text is then given to ubam_writer_header().
*/

#ifndef SAM_HEADER_H
#define SAM_HEADER_H

#include <stdbool.h>
#include <stddef.h>
#include "logger.h"
#include "strbuf.h"

/* The read groups and programs of the input files */
typedef struct {
  strbuf_t rg_lines, pg_lines;  /* @RG and @PG lines, once for each ID */
  char **rg_ids;                /* IDs of the read groups */
  size_t n_rg;
  char *last_pg;                /* ID of the last program (NULL: none) */
} sam_header_inputs_t;

/* Initialize the read groups and programs (none) */
void sam_header_inputs_init(sam_header_inputs_t *self);

/* Remove all read groups and programs */
void sam_header_inputs_clear(sam_header_inputs_t *self);

/* Free the memory */
void sam_header_inputs_free(sam_header_inputs_t *self);

/*
   Add the read groups and programs of a SAM header <text> (e.g. returned
   by bam_reader_header()), except those with an ID already added
*/
void sam_header_inputs_add(sam_header_inputs_t *self, const char *text);

/*
   Add the read groups and programs of the BAM files among the <n_files>
   files <files> (the other files are skipped)
*/
void sam_header_inputs_add_files(sam_header_inputs_t *self,
    char *const *files, size_t n_files, Logger *logger);

/*
   Append to <text> the read groups of the inputs; if <name> is not NULL
   (a destination, e.g. a sample): first a read group <name> (SM: <name>),
   then each read group X of the inputs as <name>.X, with its fields and
   SM: <name>
*/
void sam_header_append_read_groups(strbuf_t *text,
    const sam_header_inputs_t *inputs, const char *name);

/*
   Append to <text> the programs of the inputs and the @PG line of the
   program <program> (ID: <program>, or <program>.1, .2, ... if used by an
   input; PP: the last program of the inputs; VN: <version>; CL: the command
   line, see sam_header_command_line())
*/
void sam_header_append_programs(strbuf_t *text,
    const sam_header_inputs_t *inputs, const char *program,
    const char *version, const char *command_line);

/*
   The command line for the CL field (the arguments separated by spaces,
   tabs and newlines replaced by spaces); allocated, to be freed
*/
char *sam_header_command_line(int argc, char *const argv[]);

/*
   The value of the field ID of a header line (<line>, of length
   <line_len>), with its length in <len>; NULL if none
*/
const char *sam_header_line_id(const char *line, size_t line_len,
    size_t *len);

/* Unit tests (<tmpdir>: directory for temporary files) */
void sam_header_test(const char *tmpdir);

#endif
