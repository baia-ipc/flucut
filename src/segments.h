/* (c) 2026 Giorgio Gonnella */

/*
   The segments of influenza A and their recognition in the description
   of a reference sequence (the FASTA header after the name), so that
   the messages and reports can say which segment a reference is.

   The segment is taken, in this order, from:
   - "segment <n>" (n from 1 to 8, case-insensitive, e.g. NCBI records);
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

#define SEGMENTS_N 8

/* The segment (1 to 8) named in <description>, 0 if none (see above) */
int segments_of_description(const char *description);

/* The gene of segment <segment> (1 to 8), e.g. "HA" for 4 */
const char *segments_gene(int segment);

/*
   Write to <buf> (of size <size>) the label of a reference: "segment 4
   (HA), <name>" if its <description> names a segment, else <name>
*/
void segments_label(const char *name, const char *description, char *buf,
    size_t size);

/* Unit tests */
void segments_test(void);

#endif
