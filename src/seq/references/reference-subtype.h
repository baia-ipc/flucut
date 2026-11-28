/* (c) 2026 Giorgio Gonnella */

/*
   The group of a reference, parsed from its FASTA description: for
   influenza A its subtype (the H<n>N<n> of the description), for
   influenza B its lineage, where the description names it. The
   references of a group are those a sample of it is classified against
   (see auto-settings.h), and the subtype of a sample is reported from
   them.
*/

#ifndef REFERENCE_SUBTYPE_H
#define REFERENCE_SUBTYPE_H

#include <stdbool.h>
#include <stddef.h>

/*
   The subtype in a FASTA description: for influenza A the H<n>N<n>
   inside parentheses (e.g. "...(A/.../2024(H5N1)) segment 4..."),
   uppercased, or NULL where it names none; for influenza B its lineage
   (see reference-subtype.c): "Victoria", "Yamagata", or "B" where the
   description does not say which. NULL for a NULL description;
   otherwise allocated, to be freed by the caller.
*/
char* reference_parse_subtype(const char *desc);

/* The subtype named in the FASTA description <description> (the
   H<n>N<n> in it, e.g. "H5N1"), into <buf> (of size <size>); false, and
   <buf> empty, if it names none (or is NULL) */
bool reference_description_subtype(const char *description, char *buf,
    size_t size);

#endif
