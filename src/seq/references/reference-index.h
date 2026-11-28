/* (c) 2026 Giorgio Gonnella */

/*
   The references, indexed: the sequences of a FASTA file, with
   the subtype (reference-subtype.h) and the segment (segments.h)
   of each, read from its description.

   The k-mer index of the references
   ---------------------------------

   Aligning each probe to every reference costs time in proportion to the
   number of references. The references are therefore indexed by their
   k-mers (kmer-index.h of libcbaia, k = REFERENCE_KMER_LEN): the k-mers
   of a probe, and of its reverse complement, give candidate diagonals,
   and only the windows around them are aligned, with the same rules as
   the exhaustive search. k-mers occurring more than
   REFERENCE_MAX_KMER_OCC times (low complexity) are not used as seeds. A
   probe for which the k-mers give no match (e.g. an error in each k-mer,
   or a sequence absent from the references) is searched exhaustively:
   the index thus only speeds up the probes it finds. It can differ from
   the exhaustive search only where a better match than the one found
   has no exact k-mer. The exhaustive search can also be chosen for all
   probes, to check the index; or it can be left out (seeds only), for
   reads with few errors (Illumina), whose probes, if they belong to a
   reference, share exact k-mers with it, while many probes belong to no
   reference at all (e.g. reads of the host), each of which would be
   aligned to every reference.

   The index is built once and shared, read-only, by the searchers of all
   threads.
*/

#ifndef REFERENCE_INDEX_H
#define REFERENCE_INDEX_H

#include <stddef.h>
#include "fasta.h"
#include "kmer-index.h"

#define REFERENCE_KMER_LEN 10
#define REFERENCE_MAX_KMER_OCC 4096

typedef struct reference_index_t reference_index_t;

/* Index the sequences <fasta> (which must exist as long as the index is
   used) */
reference_index_t* reference_index_new(const fasta_t *fasta);

/* The sequences of the index */
const fasta_t* reference_index_fasta(const reference_index_t *self);

/* Delete (NULL is accepted) */
void reference_index_delete(reference_index_t *self);

/* The subtype of the reference <ref> (an id, -1 if its description has
   none), its number of distinct subtypes, and the name of a subtype id
   (NULL if out of range); the subtype is parsed from the FASTA
   description (e.g. the "(H5N1)" of an influenza record) */
int reference_index_subtype(const reference_index_t *self, size_t ref);
size_t reference_index_n_subtypes(const reference_index_t *self);
const char* reference_index_subtype_name(const reference_index_t *self,
    int subtype);

/* The segment of the reference <ref> (1 to 8, 0 if its description names
   none, see segments.h) */
int reference_index_segment(const reference_index_t *self, size_t ref);

/* The k-mer index of the sequences (for the searchers, reference.h) */
const kmer_index_t* reference_index_kmers(const reference_index_t *self);

#endif
