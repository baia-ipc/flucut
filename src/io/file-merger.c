/* (c) 2017-2021 Giorgio Gonnella */

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <string.h>
#include "alloc.h"
#include "file.h"
#include "heap.h"
#include "file-merger.h"

static
void get_elem_and_add_to_queue(char **elems, size_t *elem_sizes, FILE **files,
    size_t f, merge_getelem_f getelem, size_t *n_open_files, heap_t *queue,
    Logger *logger)
{
  bool has_elem;
  assert(files[f] != NULL);
  has_elem = getelem(elems + f, elem_sizes + f, files[f], logger);
  if (!has_elem)
  {
    xfclose(&files[f], logger); /* sets files[f] to NULL */
    (*n_open_files)--;
  }
  else
    heap_push(queue, elems + f);
}

static
void merge_n_files_at_once(size_t nfiles, size_t offset,
    const char *basename, const char *sfx, bool use_stdout,
    bool partial_merging, size_t partial_output_n,
    merge_getelem_f getelem, merge_writeelem_f writeelem,
    merge_cmpelems_f cmpelems, Logger *logger)
{
  char *partial_fname = NULL;
  char **filenames;
  FILE **files, *out;
  heap_t *queue;
  char **elems, **min_elem;
  size_t f, min_index, n_open_files, *elem_sizes;

  if (partial_merging) {
    partial_fname = compute_nth_filename(basename, sfx, partial_output_n);
    out = xfopen(partial_fname, NULL, "w", "merging part file", logger);
  } else {
    out = use_stdout ? stdout : xfopen(basename, sfx, "w", NULL, logger);
  }
  queue = heap_new(cmpelems, nfiles);
  XMALLOC(filenames, nfiles * sizeof(*filenames));
  XMALLOC(files, nfiles * sizeof(*files));
  XMALLOC(elems, nfiles * sizeof(*elems));
  /* zeroed: the size of an element buffer which is still NULL is 0 */
  XCALLOC(elem_sizes, nfiles, sizeof(*elem_sizes));
  n_open_files = 0;
  for (f = 0; f < nfiles; f++)
  {
    filenames[f] = compute_nth_filename(basename, sfx, f + offset);
    files[f] = xfopen(filenames[f], NULL, "r", "sorted part file", logger);
    n_open_files++;
    elems[f] = NULL;
    get_elem_and_add_to_queue(elems, elem_sizes, files, f, getelem,
        &n_open_files, queue, logger);
  }
  while (n_open_files > 0)
  {
    min_elem = heap_pop(queue);
    min_index = min_elem - elems;
    writeelem(out, *min_elem, logger);
    if (files[min_index] != NULL)
      get_elem_and_add_to_queue(elems, elem_sizes, files, min_index,
          getelem, &n_open_files, queue, logger);
  }
  heap_delete(queue);
  for (f = 0; f < nfiles; f++)
  {
    assert(files[f] == NULL);
    xfunlink(filenames[f], NULL, logger);
    free(elems[f]);
    free(filenames[f]);
  }
  free(filenames);
  free(files);
  free(elems);
  free(elem_sizes);
  if (!use_stdout)
    xfclose(&out, logger);
  if (partial_merging)
    free(partial_fname);
}

/* maximum number of files merged at once (open at the same time) */
#define NFILES_MAX 512

/* merge the files <basename><sfx>.N, for N from offset to
   offset+nfiles-1; if they are more than NFILES_MAX, they are merged in
   batches of NFILES_MAX, into the files numbered from offset+nfiles on,
   which are then merged in the same way (thus, in more than two rounds,
   if the batches are more than NFILES_MAX) */
static void merge_files_range(size_t nfiles, size_t offset,
    const char *basename, const char *sfx, bool use_stdout,
    merge_getelem_f getelem, merge_writeelem_f writeelem,
    merge_cmpelems_f cmpelems, Logger *logger)
{
  if (nfiles > NFILES_MAX)
  {
    size_t i, batch_nfiles, r_nfiles, batchnum = 0;
    for (i = 0; i < nfiles; i += NFILES_MAX) {
      r_nfiles = nfiles - i;
      batch_nfiles = (r_nfiles < NFILES_MAX) ? r_nfiles : NFILES_MAX;
      merge_n_files_at_once(batch_nfiles, offset + i, basename, sfx, false,
          true, offset + nfiles + batchnum, getelem, writeelem, cmpelems,
          logger);
      batchnum++;
    }
    merge_files_range(batchnum, offset + nfiles, basename, sfx, use_stdout,
        getelem, writeelem, cmpelems, logger);
  }
  else
  {
    merge_n_files_at_once(nfiles, offset, basename, sfx, use_stdout,
        false, 0, getelem, writeelem, cmpelems, logger);
  }
}

void merge_n_files(size_t nfiles, const char *basename, const char *sfx,
    bool use_stdout, merge_getelem_f getelem, merge_writeelem_f writeelem,
    merge_cmpelems_f cmpelems, Logger *logger)
{
  if (nfiles == 1)
  {
    if (use_stdout)
    {
      char *filename = compute_nth_filename(basename, sfx, 0);
      FILE *file = xfopen(filename, NULL, "r", NULL, logger);
      char buffer[BUFSIZ];
      size_t n;
      while ((n = fread(buffer, sizeof(char), sizeof(buffer), file)) > 0)
        xfwrite(buffer, sizeof(char), n, stdout, logger);
      LCHECK(!ferror(file), logger, "Reading file '%s' failed, %s\n",
             filename, strerror(errno));
      xfclose(&file, logger);
      xfunlink(filename, NULL, logger);
      free(filename);
    }
    else
    {
      char *filename = compute_nth_filename(basename, sfx, 0);
      xfrename(filename, NULL, basename, sfx, logger);
      free(filename);
    }
  }
  else
    merge_files_range(nfiles, 0, basename, sfx, use_stdout, getelem,
        writeelem, cmpelems, logger);
}
