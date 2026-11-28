.. _io:

File and destination I/O
------------------------

The ``io`` component contains modules for handling large files: sorting
files which do not fit in memory, managing output files for several
destinations, and computing checksums. Format-specific compression is in
:ref:`formats`.

External sorting (``file-sorter.h``, ``file-merger.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A file larger than the available memory is sorted in two steps:

1. ``sort_file_parts()`` reads the file in parts which fit in a buffer
   of a given size, splits each part into elements, sorts them, and writes
   them to the numbered temporary files ``<out><sfx>.0``, ``<out><sfx>.1``,
   ...; it returns the number of parts;
2. ``merge_n_files()`` merges the sorted parts, using a heap, into the
   output file ``<out><sfx>`` (or the standard output), deleting the
   parts; at most 512 parts are open at the same time: more parts are
   merged in batches, into further temporary files, which are then
   merged in the same way.

The elements and their order are defined by functions passed by the
caller, so that any kind of record can be sorted. For text lines, the
splitting function is ``read_lines_in_buffer()`` of ``file.h``, which
returns pointers to the lines in the buffer (each line ends with
``\n``; the lines are not 0-terminated):

.. code-block:: c

   /* compare two lines (pointers to the line pointers) */
   static int cmp_lines(const void *a, const void *b) {
     const char *x = *(char *const *)a, *y = *(char *const *)b;
     size_t lx = strcspn(x, "\n"), ly = strcspn(y, "\n");
     int c = memcmp(x, y, lx < ly ? lx : ly);
     return c != 0 ? c : (lx > ly) - (lx < ly);
   }

   static void sort_lines(char **lines, size_t n) {
     qsort(lines, n, sizeof(char *), cmp_lines);
   }

   static void output_lines(FILE *file, char **lines, size_t n,
       Logger *logger) {
     for (size_t i = 0; i < n; i++)
       xfwrite(lines[i], 1, strcspn(lines[i], "\n") + 1, file, logger);
   }

   /* merging: the elements are read with xgetline() (0-terminated, with
      the newline) */
   static bool get_line(char **line, size_t *size, FILE *file,
       Logger *logger) {
     return xgetline(line, size, file, logger);
   }

   static void write_line(FILE *file, char *line, Logger *logger) {
     (void)logger;
     fputs(line, file);
   }

   static int cmp_line_ptrs(const void *a, const void *b) {
     return strcmp(*(char *const *)a, *(char *const *)b);
   }

   /* sort big.txt into big.txt.sorted, using a buffer of 1 GiB */
   size_t n_parts = sort_file_parts("big.txt", NULL, "big.txt", ".sorted",
       (size_t)1 << 30, 1024, read_lines_in_buffer, sort_lines,
       output_lines, logger);
   merge_n_files(n_parts, "big.txt", ".sorted", false, get_line,
       write_line, cmp_line_ptrs, logger);

The buffer must be larger than the longest line (otherwise the program
exits with a message suggesting a larger buffer).

Checksums (``md5.h``)
~~~~~~~~~~~~~~~~~~~~~

The MD5 digest (RFC 1321) of data or of a file, as 32 lowercase
hexadecimal digits, e.g. for the checksums of files in metadata (such
as the ``md5`` of the onlists of seqspec files) or for checking
downloads; not for security purposes:

.. code-block:: c

   char hex[MD5_HEX_LEN + 1];
   md5_file("barcodes.txt.gz", hex, logger);   /* as md5sum */

   md5_t md5;                                  /* data in parts */
   md5_init(&md5);
   md5_update(&md5, data, len);
   md5_update(&md5, more, more_len);
   md5_final(&md5, hex);

.. _outfiles:

Output files of several destinations (``outfiles.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A demultiplexer, or any program distributing the reads to several
*destinations* (e.g. the samples), writes a file for each of them. An
``outfiles_t`` manages these files:

* the destinations are numbered from 0; each has a name, used in the
  file name;
* for an input file, the files are
  ``<outdir>/<input basename>.<destination>.fastq.gz`` (the suffix
  of a file of reads of the input, ``.fastq``, ``.fq`` (optionally
  ``.gz``) or ``.bam``, is removed); for a
  group of input files (e.g. the files of a barcode directory),
  ``<outdir>/<group>/<destination>.fastq.gz``;
* a file is created when the first data are written to it, so that
  destinations without reads have no file;
* the data are written already compressed (gzip members, e.g. by a
  ``gzip_compressor_t``), so that they can be compressed in parallel;
* with ``outfiles_set_append()``, the data are appended to existing
  files (e.g. when the files of a run are processed one at a time);
* BAM output: with ``outfiles_set_bam()``, the files are ``.bam``
  instead of ``.fastq.gz``; the data are BGZF blocks; the header (given,
  BGZF-compressed, for each destination) is written to each new file,
  and the end-of-file marker when it is closed (and removed before
  appending to an existing file).

.. code-block:: c

   outfiles_t *out = outfiles_new("out", "run1.fastq.gz", n_samples + 1,
                                  logger);
   ...
   /* compressed: a gzip member with the records of destination d */
   outfiles_write(out, d, names[d], compressed.s, compressed.len);
   ...
   outfiles_delete(out);   /* closes the files */

The output directory (and the group directory) is created if needed.
:ref:`recipe-demux` shows a complete program.
