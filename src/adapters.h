/* (c) 2026 Giorgio Gonnella */

/*
   The sequencing adapters of Oxford Nanopore, searched at the ends of a
   read at which no primer was found (see amplicon.h), so that the
   adapter, the barcode and what lies outside them can be removed.

   A read of a ligation library begins with the adapter (with native
   barcoding: the ligation adapter, then the barcode flanked by two
   fixed sequences) and, if it was sequenced to the end, ends with the
   reverse complement of the adapter of the other end of the molecule.
   When the primer of an end is found, everything outside it is removed
   with it; when it is not (the read is broken, or begins too far inside
   the primer), the adapter is what can still be recognised.

   The kits (ADAPTERS_KITS), with their sequences from Dorado
   (github.com/nanoporetech/dorado, commit 8b8fc5d:
   dorado/utils/barcode_kits.cpp and dorado/demux/adapter_primer_kits.cpp):

   - nb: native barcoding (e.g. SQK-NBD114-24/96, EXP-NBD114): at the
     start of the read, ATTGCTAAGGTTAA, the barcode (24 bases),
     CAGCACCT; at the end, the reverse complement of ATTGCTAAGGTTAA,
     the barcode, CAGCACC;
   - rbk: rapid barcoding (e.g. SQK-RBK114-24/96): at the start of the
     read only, GCTTGGGTGTTTAACC, the barcode (24 bases), and the 50
     bases of the rapid adapter;
   - lsk: ligation without barcodes (e.g. SQK-LSK114): at the start of
     the read CCTGTACTTCGTTCAGTTACGTATTGC, at the end
     AGCAATACGTAACTGAAC.

   The barcode is searched as Ns, which match any base: any barcode of
   the kit is recognised, and none is identified (demultiplexing is not
   the purpose: the reads normally come already demultiplexed). The
   errors allowed are ADAPTERS_ERROR_RATE of the fixed bases (not the
   barcode), and an adapter without a barcode may overhang the read by
   up to half of its length (the read begins or ends inside it, see
   adaptor-search.h; with a barcode, the part left in the read could be
   mostly the barcode, which matches any sequence: an adapter with a
   barcode must be in the read, its outer bases missing counted as
   errors).

   An adapter is accepted only if it begins (start of the read) or ends
   (end of the read) within <window> bases of that end of the read, as
   the primers (amplicon.h); of the adapters of the kits found at an end
   the innermost is taken (e.g. with native barcoding, the barcode
   adapter, not the ligation adapter before it).
*/

#ifndef ADAPTERS_H
#define ADAPTERS_H

#include <stdbool.h>
#include <stddef.h>
#include "logger.h"

/* the kits, for the option --adapters (see above) */
#define ADAPTERS_KITS "nb, rbk, lsk"

/* the default: all the kits */
#define ADAPTERS_DEFAULT "nb,rbk,lsk"

/* the errors allowed, as a fraction of the fixed bases of an adapter */
#define ADAPTERS_ERROR_RATE 0.2

/* the length of the barcodes */
#define ADAPTERS_BARCODE_LEN 24

/* an adapter found at an end of a read */
typedef struct {
  bool found;
  size_t kit;        /* see adapters_kit_name() */
  long start, end;   /* 0-based, inclusive, in the read */
  long n_errors;
} adapter_hit_t;

typedef struct adapters_t adapters_t;

/*
   The adapters of the kits <kits> (a comma-separated list of names, see
   above; "none": no kit, see adapters_any()); the program is exited with
   an error for an unknown kit. <window>: see above. The object contains
   buffers: use one for each thread.
*/
adapters_t* adapters_new(const char *kits, long window, Logger *logger);

/* Whether the list of kits <kits> is valid (see adapters_new()) */
bool adapters_valid(const char *kits);

/* Whether any kit is searched */
bool adapters_any(const adapters_t *self);

/* Search the adapters at the start of the read <seq> (length <len>) */
void adapters_search_start(adapters_t *self, const char *seq, size_t len,
    adapter_hit_t *hit);

/* Search the adapters at the end of the read <seq> (length <len>) */
void adapters_search_end(adapters_t *self, const char *seq, size_t len,
    adapter_hit_t *hit);

/* The name of the kit number <kit> (e.g. "nb") */
const char* adapters_kit_name(size_t kit);

/* Delete (NULL is accepted) */
void adapters_delete(adapters_t *self);

/* Unit tests */
void adapters_test(void);

#endif
