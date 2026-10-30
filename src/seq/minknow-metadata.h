/* (c) 2026 Giorgio Gonnella */

/*
   Read metadata of MinKNOW: the FASTQ files written by MinKNOW contain the
   read metadata in the title as fields key=value, separated by spaces
   (e.g. runid=... ch=12 start_time=...), while Dorado writes SAM tags
   (see sam-tags.h). The fields can be converted to SAM tags:
   runid: ri:Z, read: rn:i, ch: ch:i, start_time: st:Z, flow_cell_id:
   PU:Z, protocol_group_id: pg:Z, sample_id: si:Z, barcode: ob:Z,
   barcode_alias: oa:Z, parent_read_id: pi:Z, basecall_model_version_id:
   bm:Z.

   rn, ch, st and pi are the tags of Dorado; PU (platform unit, a standard
   tag) contains the flow cell, as in the read groups of Dorado; the other
   tags are not standard (the SAM specification reserves the tags with a
   lowercase letter for end users) and avoid the tags written by Dorado
   and minimap2. The integer fields are written as Z if their value is not
   an integer.
*/

#ifndef MINKNOW_METADATA_H
#define MINKNOW_METADATA_H

#include <stdbool.h>
#include <stddef.h>
#include "strbuf.h"

/* SAM tag of the ONT barcode (barcode= of MinKNOW; Dorado writes it in
   the BC tag, which is the standard tag of the sample barcode) */
#define MINKNOW_BARCODE_SAM_TAG "ob"

/*
   If <field> (of length <len>) is a field key=value of the read metadata
   of MinKNOW (see above), append it to <out> as a SAM tag,
   preceded by a tab (e.g. "ch=12" as "\tch:i:12"), and return true;
   otherwise return false (and <out> is unchanged)
*/
bool minknow_metadata_append_sam_tag(strbuf_t *out, const char *field,
    size_t len);

/* Unit tests */
void minknow_metadata_test(void);

#endif
