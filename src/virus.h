/* (c) 2026 Giorgio Gonnella */

/*
   The type of the influenza virus whose genome was amplified, A or B.

   The two types are amplified by the same multisegment RT-PCR, with
   different universal primers, because their segment termini differ
   (see primers.h). They also number their segments differently: the
   largest segment of influenza A is PB2, that of influenza B is PB1
   (see segments.h). Everything else flucut does is the same for both.

   A run can hold samples of either type, so the type is not a property
   of the run: it is a property of a reference (read from its
   description) and of a sample (the type of the references its reads
   were located in). What the user chooses with --virus is which
   primers are searched and which references are used (see options.h):
   both types (the default), or one of them alone.
*/

#ifndef VIRUS_H
#define VIRUS_H

#include <stdbool.h>
#include <stddef.h>

typedef enum {
  VIRUS_A,
  VIRUS_B,
} virus_type_t;

#define VIRUS_N_TYPES 2

/* The name of a type: "A" or "B" */
const char* virus_type_name(virus_type_t type);

/*
   The type given by the user as <s>: "A" or "B", in any case, also
   written "flu A", "influenza B" or "FluB"; false, and <type>
   unchanged, if <s> names neither.
*/
bool virus_type_parse(const char *s, virus_type_t *type);

/*
   Whether <description> (the FASTA description of a reference, or the
   name of an IRMA seed) names influenza B: "influenza B" or
   "influenzavirus B" in any case, a strain name beginning with "B/"
   (e.g. "B/Brisbane/60/2008"), or a name beginning with "B_" (the IRMA
   seeds B_HA, B_PB1, ...). A description which names neither type is
   taken for influenza A, as the records of influenza A outnumber
   those of influenza B and the A records are the ones which may omit
   the type (e.g. "A/Puerto Rico/8/1934" is also recognised, but
   "segment 4 hemagglutinin" alone gives A).
*/
bool virus_description_is_b(const char *description);

/* The type named in <description>, VIRUS_A if none is (see above) */
virus_type_t virus_type_of_description(const char *description);

/*
   A set of types: which types a run may hold (--virus). The primers of
   the types in the set are searched and the references of those types
   are used; the type of each sample is then read from the references
   its reads were located in.
*/
typedef unsigned virus_set_t;

#define VIRUS_SET_OF(TYPE) ((virus_set_t)1 << (TYPE))
#define VIRUS_SET_A VIRUS_SET_OF(VIRUS_A)
#define VIRUS_SET_B VIRUS_SET_OF(VIRUS_B)
#define VIRUS_SET_BOTH (VIRUS_SET_A | VIRUS_SET_B)

/* Whether <set> holds <type> */
bool virus_set_has(virus_set_t set, virus_type_t type);

/*
   The set given by the user as <s>: a type ("A", "B", see
   virus_type_parse()) or "both" (also "auto" and "ab"); false, and
   <set> unchanged, if <s> names none.
*/
bool virus_set_parse(const char *s, virus_set_t *set);

/* The name of a set, for the messages: "A", "B" or "A and B" */
const char* virus_set_name(virus_set_t set);

/*
   The type of a set which holds one type, and VIRUS_A for the set of
   both: the type whose numbering of the segments is used where the
   user gives a segment by its number and no reference says which type
   is meant (see segments.h).
*/
virus_type_t virus_set_main(virus_set_t set);

/* Unit tests */
void virus_test(void);

#endif
