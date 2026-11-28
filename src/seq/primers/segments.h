/* (c) 2026 Giorgio Gonnella */

/*
   The segments of influenza and their recognition in the description
   of a reference sequence (the FASTA header after the name), so that
   the messages and reports can say which segment a reference is.

   A segment is held as its gene: 1 to 8 are PB2, PB1, PA, HA, NP, NA,
   M and NS, which is the order in which influenza A numbers them and
   the order in which the reports show them, for both types. Influenza
   B numbers its two largest segments the other way round (its segment
   1 is PB1 and its segment 2 is PB2), so the number a description
   gives and the number shown to the user are translated by
   segments_number() and segments_of_number(); every other segment has
   the same number in both types.

   The segment is taken, in this order, from:
   - "segment <n>" (n from 1 to 8, case-insensitive, e.g. NCBI
     records), read in the numbering of the type the description names
     (see virus.h);
   - the first word which names a gene or its protein: PB2, PB1, PA, HA
     or hemagglutinin, NP or nucleoprotein (nucleocapsid), NA or
     neuraminidase, M, M1, M2, MP (IRMA) or matrix, NS, NS1, NS2, NEP or
     nonstructural (non-structural), the gene names in upper case (so
     that e.g. "na" or "pa" in a text is not taken for one), the
     protein names in any case; words are separated by any
     character which is not a letter or a digit, so that e.g. the IRMA
     name A_HA_H5 is recognised.
   A description naming none of these gives no segment.
*/

#ifndef SEGMENTS_H
#define SEGMENTS_H

#include <stddef.h>
#include "virus.h"

/* the number of segments (SEGMENTS_N) and each of them by its gene
   (SEGMENT_PB2 to SEGMENT_NS, numbered as in influenza A), written by
   the build from metadata/segments.yaml, which also lists the words
   which name a gene in the description of a reference */
#include "segments-table.h"

/* The segment (1 to 8) named in <description>, 0 if none (see above) */
int segments_of_description(const char *description);

/* The gene of segment <segment> (1 to 8), e.g. "HA" for 4 */
const char *segments_gene(int segment);

/* The number by which <type> names segment <segment> (1 to 8), and the
   segment of the number <number> in <type>: the identity for influenza
   A, and PB2 and PB1 exchanged for influenza B (see above) */
int segments_number(int segment, virus_type_t type);
int segments_of_number(int number, virus_type_t type);

/* The segment (1 to 8) given by the user as <s>: its number ("4", in
   the numbering of <type>) or a name of its gene or protein, in any
   case ("HA", "ha", "hemagglutinin"); 0 if <s> is neither */
int segments_parse(const char *s, virus_type_t type);

/* Write to <buf> (of size <size>) the name of segment <segment> (1 to
   8) in the numbering of <type>, e.g. "segment 4 (HA)" */
void segments_name(int segment, virus_type_t type, char *buf,
    size_t size);

/*
   Write to <buf> (of size <size>) the name of segment <segment> in the
   numbering of the type named in <description> (see virus.h), e.g.
   "segment 1 (PB1)" for a reference of influenza B. <description> may
   be NULL (the numbering of influenza A is then used).
*/
void segments_name_in(int segment, const char *description, char *buf,
    size_t size);

/*
   Write to <buf> (of size <size>) the label of a reference: "segment 4
   (HA), <name>" if its <description> names a segment, else <name>
*/
void segments_label(const char *name, const char *description, char *buf,
    size_t size);

/* Unit tests */
void segments_test(void);

#endif
