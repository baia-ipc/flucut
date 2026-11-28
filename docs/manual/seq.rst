.. _seq:

Sequences and sequencing reads
------------------------------

The ``seq`` component contains sequence algorithms and the
format-independent ``reads-reader`` interface. The modules for file
formats are described in :ref:`formats`, MinKNOW metadata in
:ref:`nanopore`, and single-cell barcode and complexity modules in
:ref:`singlecell`.

Sequence utilities (``sequtil.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``void reverse_complement(const char *seq, size_t len, char *dest)``
   write the reverse complement of a sequence to ``dest`` (at least
   ``len + 1`` bytes, 0-terminated); IUPAC codes are complemented, ``U``
   is complemented to ``A``, the case is kept and other characters are
   unchanged;
``void py_slice(long len, long start, long stop, long *from, long *slice_len)``
   the bounds of the Python slice ``s[start:stop]`` of a string of
   length ``len`` (negative values count from the end, values out of the
   range are clipped), useful for reproducing exactly the results of
   Python programs.

K-mer index (``kmer-index.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The k-mer index finds where a query may align in a set of sequences
(e.g. reference sequences), without aligning it to the whole of each:
the k-mers of the query which occur in the indexed sequences give
**candidate diagonals** (position in the sequence minus position in the
query), around which the caller aligns the query (e.g. with edlib).

The k-mers are encoded with 2 bits per base (k up to
``KMER_INDEX_MAX_K``, 12); k-mers with a character other than A, C, G, T
are not indexed. The table of the offsets has 4^k + 1 entries (4 MB for
k = 10). The index is only read after its construction, and can be
shared by threads; the buffer of the candidates belongs to the caller
(one per thread). A k-mer occurring more often than a limit (e.g. of low
complexity) can be skipped. A candidate needs an exact k-mer: a query
with an error in each of its k-mers gives none.

.. code-block:: c

   kmer_index_t *ix = kmer_index_new(seqs, lens, n_seqs, 10);
   kmer_index_hit_t *hits = NULL;
   size_t cap = 0;
   size_t n = kmer_index_diagonals(ix, query, query_len, 1000, &hits,
                                   &cap);
   for (size_t i = 0; i < n; i++)
     align(query, seqs[hits[i].seq], hits[i].diagonal);
   free(hits);
   kmer_index_delete(ix);

The candidates are sorted by sequence and diagonal, without duplicates,
so that nearby diagonals can be merged into one window to align.

FASTQ or BAM (``reads-reader.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most programs accept both formats. ``reads_reader_t`` opens either,
recognized from the content, with the interface of the FASTQ reader:

.. code-block:: c

   reads_reader_t *reader = reads_reader_new(filename, logger);
   while ((n = reads_reader_read_batch(reader, &batch, 1000)) > 0)
     ...
   reads_reader_delete(reader);

Further functions:

``bool reads_reader_is_bam(const reads_reader_t *self)``
   whether the file read is a BAM file (e.g. to read its header with a
   ``bam_reader_t``, or to write BAM output);
``bool reads_reader_is_reads_filename(const char *name)``
   whether a name is that of a file of reads: suffix ``.fastq``, ``.fq``
   (optionally ``.gz``) or ``.bam``; e.g. as a filter of ``list_dir()``;
``size_t reads_reader_filename_suffix_len(const char *name)``
   the length of that suffix at the end of the name (0 if none), e.g. to
   remove it;
``size_t reads_reader_count(const char *filename, Logger *logger)``
   the number of reads of a file (e.g. for the total of a progress bar);
   the file is read and parsed (thus also multi-line FASTQ records are
   counted correctly), which takes about as long as decompressing it; the
   program is exited if the file is not valid;
``bool reads_reader_file_complete(const char *path)``
   whether a file written by another program is complete: a gzipped file
   is not truncated and ends with a newline, a FASTQ file ends with a
   newline, a BAM file ends with the BGZF end-of-file marker. A gzipped
   file is decompressed completely (its last character is not known
   otherwise), which takes the time of reading it; an uncompressed or BAM
   file is read only at its start and end. This is used by programs
   which process the files of a sequencing run while they are written
   (e.g. in real-time sequencing workflows), together with the age of the
   file:

   .. code-block:: c

      struct stat st;
      if (stat(path, &st) == 0 &&
          difftime(time(NULL), st.st_mtime) >= min_age &&
          reads_reader_file_complete(path))
        process(path);
