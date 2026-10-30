/* (c) 2026 Giorgio Gonnella */

/*
   SAM tags in text format, in read titles: the title of a read (FASTQ
   header without '@', or the title returned by the BAM reader) is the
   read name followed by fields; in the FASTQ written by Dorado, and in the
   output of programs which keep the tags (e.g. for minimap2 -y or samtools
   import -T '*'), the fields are SAM tags in text format (e.g. ch:i:12),
   separated by tabs. The uBAM writer encodes them in binary form. For the
   read metadata of MinKNOW (fields key=value), see minknow-metadata.h.
*/

#ifndef SAM_TAGS_H
#define SAM_TAGS_H

#include <stdbool.h>
#include <stddef.h>
#include "strbuf.h"

/*
   Is the text <field> (of length <len>) a tag in SAM text format, i.e.
   a two-character name (a letter and a letter or digit), a type (A, i, f,
   Z, H or B) and the value, separated by ':' (e.g. "RG:Z:run1")
*/
bool sam_tag_is_valid(const char *field, size_t len);

/*
   Append to <out> a tab and the start of a tag <tag> (two characters) of
   type <type> (e.g. "\tCB:Z:"), to which the value is then appended
*/
void sam_tag_start(strbuf_t *out, const char *tag, char type);

/* Append to <out> a tab and the string tag <tag>:Z:<value> (<len> bytes) */
void sam_tag_append_string(strbuf_t *out, const char *tag,
    const char *value, size_t len);

/* Append to <out> a tab and the integer tag <tag>:i:<value> */
void sam_tag_append_int(strbuf_t *out, const char *tag, long value);

/* Unit tests */
void sam_tags_test(void);

#endif
