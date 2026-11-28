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

FASTA files (``fasta.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~

``fasta_read()`` reads a plain FASTA file into memory. Each record has a
name (the first word of the header), the rest of the header as its
description, and an upper-case sequence. IUPAC codes are accepted; an
empty file, an empty record, or a sequence with another character is
reported as an error. This is intended for small sequence collections,
such as reference sequences.

IUPAC codes (``iupac.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~

``iupac_bases(code)`` returns the canonical bases represented by an
IUPAC code, and ``iupac_is_seq(seq, len)`` checks that a non-empty
sequence contains only IUPAC codes. ``iupac_equal_chars()`` provides the
characters considered equivalent for matching, including ambiguity
codes and their bases. The IUPAC pair type is defined by libcbaia and
does not depend on an alignment library.

Approximate pattern scan (``pattern-scan.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The pattern scanner finds approximate matches of a short IUPAC pattern
in one pass over a sequence, using a Myers bit-parallel scan for patterns
up to one machine word and a linear-space fallback for longer patterns.
It can report every match or the best forward and reverse matches.

.. code-block:: c

   pattern_scan_t *scan = pattern_scan_new(pattern, pattern_len);
   pattern_scan_end_t ends[PATTERN_SCAN_MAX_MATCHES];
   size_t n = pattern_scan_run(scan, seq, seq_len, max_errors, ends,
                               PATTERN_SCAN_MAX_MATCHES);
   pattern_scan_delete(scan);

Read statistics (``read-stats.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``read_stats_t`` keeps read counts, bases, length histograms, and quality
statistics in constant memory. Initialize it with the number of sets and
their names, then add reads to a set by index. The set names remain owned
by the caller and must outlive the statistics object. Summaries include
mean, median, and N50 lengths and the median Phred score of the mean base
error probability.

.. code-block:: c

   static const char *const names[] = {"before", "after"};
   read_stats_t stats;
   read_stats_init(&stats, 2, names);
   read_stats_add(&stats, 0, len, qualities);
   read_stats_summary_t summary;
   read_stats_summarize(&stats, 0, &summary);

Read-title metadata (``read-metadata.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``read_metadata_add_title()`` extracts run, flow-cell, model, sample,
experiment, instrument, and start-time fields from Oxford Nanopore and
Illumina read titles. Values that differ across reads are reported as
``several``; the earliest start time is kept. The module can merge
metadata collected from several batches and save or reload it as text.

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
