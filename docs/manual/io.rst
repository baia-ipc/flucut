.. _io:

Large files and compression
---------------------------

The ``io`` component contains modules for handling large files: sorting
files which do not fit in memory, compressing output data in parallel,
and computing checksums.

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

Parallel compression (``gzip-compressor.h``, ``bgzf-compressor.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compressing the output is often the most expensive part of a read
processing program. The two compressors make it parallel, using a
property of the gzip format: **the concatenation of gzip members is a
valid gzip file**. Each block of output (e.g. the records of a batch) is
compressed independently, as a gzip member, by the thread which produced
it; the compressed blocks are then written to the file in the right
order (see :ref:`recipe-threads`).

.. code-block:: c

   gzip_compressor_t *gz = gzip_compressor_new(6);   /* level 1-9 */
   strbuf_t out;
   strbuf_init(&out);
   gzip_compressor_run(gz, text.s, text.len, &out);  /* appended to out */
   fwrite(out.s, 1, out.len, file);                  /* a gzip member */
   strbuf_free(&out);
   gzip_compressor_delete(gz);

Each thread uses its own compressor (the zlib stream state is reused
for all its blocks).

``bgzf-compressor.h`` does the same in the BGZF format of BAM files:
blocks of at most 64 KiB, each a gzip member with an extra field giving
its compressed size (so that BAM readers can seek). A BGZF file ends
with an empty block, the end-of-file marker ``bgzf_eof``
(``BGZF_EOF_LEN`` bytes), written after the last block:

.. code-block:: c

   bgzf_compressor_t *bgzf = bgzf_compressor_new(6);
   bgzf_compressor_run(bgzf, bam_data.s, bam_data.len, &out);
   ...
   fwrite(bgzf_eof, 1, BGZF_EOF_LEN, file);

The output files of several destinations (``outfiles.h``,
:ref:`outfiles`) write such compressed blocks, and add the BGZF
end-of-file marker to BAM files.

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
