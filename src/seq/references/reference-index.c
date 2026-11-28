/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "fasta.h"
#include "kmer-index.h"
#include "reference-index.h"
#include "reference-subtype.h"
#include "segments.h"

struct reference_index_t {
  const fasta_t *fasta;
  kmer_index_t *kmers;
  /* the subtype of each reference (an id into <subtype_names>, -1 if its
     description has none), parsed from the FASTA description (e.g. the
     "(H5N1)" of an influenza record): the references of different subtypes
     can be told apart, so that a read is classified against the subtype it
     belongs to (see reference_set_restrict) */
  int *subtype;
  char **subtype_names;
  size_t n_subtypes;
  /* the segment of each reference (1 to 8, 0 if its description names
     none, see segments.h), for reference_set_restrict_segment() */
  int *segment;
};

/* fill self->subtype (an id per record) and self->subtype_names */
static void index_subtypes(reference_index_t *self) {
  size_t n = self->fasta->n_records;
  XMALLOC(self->subtype, sizeof(int) * (n > 0 ? n : 1));
  self->subtype_names = NULL;
  self->n_subtypes = 0;
  size_t cap = 0;
  for (size_t i = 0; i < n; i++) {
    char *name =
        reference_parse_subtype(self->fasta->records[i].description);
    if (name == NULL) {
      self->subtype[i] = -1;
      continue;
    }
    int id = -1;
    for (size_t j = 0; j < self->n_subtypes; j++)
      if (strcmp(self->subtype_names[j], name) == 0) {
        id = (int)j;
        break;
      }
    if (id < 0) {
      if (self->n_subtypes == cap) {
        cap = cap == 0 ? 8 : cap * 2;
        XREALLOC(self->subtype_names, sizeof(char *) * cap);
      }
      self->subtype_names[self->n_subtypes] = name;
      id = (int)self->n_subtypes++;
    } else {
      free(name);
    }
    self->subtype[i] = id;
  }
}

reference_index_t* reference_index_new(const fasta_t *fasta) {
  assert(fasta != NULL);
  reference_index_t *self;
  XMALLOC(self, sizeof(reference_index_t));
  self->fasta = fasta;
  index_subtypes(self);
  const char **seqs;
  size_t *lens;
  size_t n = fasta->n_records;
  XMALLOC(self->segment, sizeof(int) * (n > 0 ? n : 1));
  for (size_t i = 0; i < n; i++)
    self->segment[i] = fasta->records[i].description == NULL ? 0 :
        segments_of_description(fasta->records[i].description);
  XMALLOC(seqs, sizeof(char *) * (n > 0 ? n : 1));
  XMALLOC(lens, sizeof(size_t) * (n > 0 ? n : 1));
  for (size_t i = 0; i < n; i++) {
    seqs[i] = fasta->records[i].seq;
    lens[i] = fasta->records[i].seq_len;
  }
  self->kmers = kmer_index_new(seqs, lens, n, REFERENCE_KMER_LEN);
  free(seqs);
  free(lens);
  return self;
}

const fasta_t* reference_index_fasta(const reference_index_t *self) {
  return self->fasta;
}

void reference_index_delete(reference_index_t *self) {
  if (self == NULL)
    return;
  kmer_index_delete(self->kmers);
  for (size_t i = 0; i < self->n_subtypes; i++)
    free(self->subtype_names[i]);
  free(self->subtype_names);
  free(self->subtype);
  free(self->segment);
  free(self);
}

int reference_index_subtype(const reference_index_t *self, size_t ref) {
  return self->subtype[ref];
}

size_t reference_index_n_subtypes(const reference_index_t *self) {
  return self->n_subtypes;
}

int reference_index_segment(const reference_index_t *self, size_t ref) {
  return self->segment[ref];
}

const char* reference_index_subtype_name(const reference_index_t *self,
    int subtype) {
  if (subtype < 0 || (size_t)subtype >= self->n_subtypes)
    return NULL;
  return self->subtype_names[subtype];
}

const kmer_index_t* reference_index_kmers(const reference_index_t *self) {
  return self->kmers;
}
