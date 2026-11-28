.. _formats:

File formats
------------

The ``formats`` tree groups modules by the file format they read or
write. Each format has its own source directory:

* ``formats/gzip``: gzip compression;
* ``formats/bgzf``: BGZF compression used by BAM;
* ``formats/json``: JSON writing;
* ``formats/sam_bam``: SAM headers and tags, BAM reading, and unaligned
  BAM writing;
* ``formats/fastq``: FASTQ reading and writing.

The format-specific APIs and examples follow. The format-independent
``reads-reader`` interface is described in :ref:`seq`; the generic
destination manager ``outfiles`` is in :ref:`io`.

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


Reading FASTQ files (``fastq-reader.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The FASTQ reader reads plain or gzipped files (recognized from the
content) in **batches** of records. The parsing follows the Biopython
``FastqGeneralIterator``: the title is kept as it is (except trailing
whitespace), sequences and qualities can span several lines, and the
title of the ``+`` line is optional, but must be equal to that of the
record, if present. Unlike Biopython, which skips any text before the
first record, a file whose first line does not start with ``@`` is not
valid. Invalid records interrupt the program with a message giving the
file, the line and the problem.

A batch (``fastq_batch_t``) stores all its records in one buffer
(``data``); each record (``fastq_record_t``) gives the offsets and
lengths of its parts in the buffer:

.. code-block:: c

   fastq_reader_t *reader = fastq_reader_new("reads.fastq.gz", logger);
   fastq_batch_t batch;
   fastq_batch_init(&batch);
   size_t n;
   while ((n = fastq_reader_read_batch(reader, &batch, 1000)) > 0) {
     for (size_t i = 0; i < n; i++) {
       const fastq_record_t *r = &batch.records[i];
       const char *title = batch.data.s + r->title_from;  /* without @ */
       const char *seq = batch.data.s + r->seq_from;
       const char *qual = batch.data.s + r->qual_from;    /* seq_len chars */
       process(title, r->title_len, seq, qual, r->seq_len);
     }
   }
   fastq_batch_free(&batch);
   fastq_reader_delete(reader);

The strings are not 0-terminated: use the lengths. A batch replaces its
content at each call, and keeps its memory, so that reading a file
allocates memory only for the first batches. In a multithreaded program,
each thread reads its batch (under a lock) and processes it, see
:ref:`recipe-threads`. ``fastq_reader_count()`` returns the number of
records of a file (reading it as above, without the log message).

Reading BAM files (``bam-reader.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The BAM reader reads BAM files (e.g. the unaligned BAM of Dorado and
MinKNOW) without htslib, into the same batches of FASTQ-like records:

* the title is the read name, followed by the tags in SAM text format,
  separated by tabs (as in the FASTQ output of Dorado), e.g.
  ``read1<TAB>RG:Z:run1<TAB>ch:i:42``; array tags (type ``B``, e.g. the
  move table) and the base modification tags (``MM``, ``ML``, ``MN``) are
  not included;
* the sequence and the qualities are those of the original read (reverse
  complemented if the record is mapped to the minus strand); quality
  values above 93 are written as ``~`` (the maximum of SAM); missing
  qualities are written as ``!`` and the record is marked
  (``qual_missing``);
* secondary and supplementary alignments are skipped.

``bam_reader_header()`` returns the SAM header text (e.g. to copy the
``@RG`` lines to the output); ``bam_reader_is_bam()`` tells whether a
file is a BAM file (by its content); ``bam_reader_count()`` counts the
records; ``bam_reader_has_eof_marker()`` tells whether the file ends with
the BGZF end-of-file marker, i.e. it was completely written.

Writing FASTQ and BAM records (``fastq-writer.h``, ``ubam-writer.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The writers append records to a ``strbuf_t``, which is then compressed
(with the gzip or BGZF compressor) and written to a file:

``fastq_writer_record(&out, title, title_len, seq, qual, len)``
   a FASTQ record, in four lines (the title without ``@``);
``ubam_writer_header(&out, text)``
   the header of an unaligned BAM file, with the SAM header text (no
   reference sequences; see :ref:`sam-header`);
``ubam_writer_record(&out, title, title_len, seq, qual, len)``
   an unmapped BAM record (flag 4): the title contains the read name and
   SAM tags in text format, separated by tabs (as returned by the BAM
   reader); the tags are encoded in binary form (integers with the
   smallest type, as htslib does); fields which are not valid SAM tags
   are not written, nor repeated tags (a BAM record cannot have two tags
   with the same name: the first valid one is kept). The read name is the first word of the title
   (leading whitespace is skipped; ``*``, i.e. ``UBAM_NO_READ_NAME``, if
   the title is empty); names longer than ``UBAM_MAX_READ_NAME`` (254)
   are truncated. Quality characters below ``!`` or above ``~`` (out of
   the range of SAM) are clamped to it. Missing qualities (``0xFF`` in
   BAM) are written if ``qual`` is ``NULL``; the BAM reader gives them as
   ``!`` (as FASTQ has no missing qualities) and sets ``qual_missing`` of
   the record, so that a record read and written again (with ``NULL``
   qualities if ``qual_missing``) keeps them missing.

Thus, a program can read FASTQ or BAM, and write FASTQ or BAM, keeping
the tags of the input (e.g. the read group, the channel, or the tags
added by the program, such as ``CB`` and ``UB`` for single-cell data).

.. _sam-header:

The header of BAM output (``sam-header.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A program which writes BAM files from BAM input (e.g. the unaligned BAM
of Dorado) keeps the provenance of the reads in the SAM header text given
to ``ubam_writer_header()``: the read groups (``@RG``) of the inputs, so
that the ``RG`` tags of the records still refer to a read group of the
header, and the programs (``@PG``) of the inputs, followed by the
``@PG`` line of the program itself, linked to the last one of the inputs
(``PP``), as samtools does.

The read groups and programs of the inputs are collected in a
``sam_header_inputs_t`` (a plain structure: ``rg_ids`` and ``n_rg`` are
the IDs of the read groups, ``last_pg`` the ID of the last program, or
``NULL``):

``void sam_header_inputs_init(sam_header_inputs_t *self)``, ``void sam_header_inputs_clear(sam_header_inputs_t *self)``, ``void sam_header_inputs_free(sam_header_inputs_t *self)``
   initialize it (no read groups and programs), remove all read groups
   and programs, free its memory;
``void sam_header_inputs_add(sam_header_inputs_t *self, const char *text)``
   add the ``@RG`` and ``@PG`` lines of a SAM header text (e.g. returned
   by ``bam_reader_header()``), except those with an ID already added
   (e.g. the same read group in several files of a run);
``void sam_header_inputs_add_files(sam_header_inputs_t *self, char *const *files, size_t n_files, Logger *logger)``
   the same, for the BAM files among the input files (the other files,
   e.g. FASTQ, are skipped).

The header text is then built with:

``void sam_header_append_read_groups(strbuf_t *text, const sam_header_inputs_t *inputs, const char *name)``
   append the read groups of the inputs; with a ``name`` (a destination,
   e.g. a sample, of a program which distributes the reads), first a read
   group ``name`` (``SM:name``), then each read group ``X`` of the inputs
   renamed as ``name.X``, with its fields and ``SM:name``; the program
   then writes ``RG:Z:name.X`` in the records of the destination (or
   ``RG:Z:name``, for reads without read group);
``void sam_header_append_programs(strbuf_t *text, const sam_header_inputs_t *inputs, const char *program, const char *version, const char *command_line)``
   append the programs of the inputs and the ``@PG`` line of the
   program: ``ID:program`` (``program.1``, ``.2``, ... if an input was
   already written by it), ``PN:program``, ``PP:`` the last program of the
   inputs, ``VN:version``, ``CL:command_line``;
``char *sam_header_command_line(int argc, char *const argv[])``
   the command line for the ``CL`` field: the arguments separated by
   spaces, with tabs and newlines replaced by spaces (allocated);
``const char *sam_header_line_id(const char *line, size_t line_len, size_t *len)``
   the value of the ``ID`` field of a header line (not 0-terminated, its
   length in ``*len``); ``NULL`` if none.

For example, a demultiplexer writing a BAM file for each sample sets the
header of each destination of the output files (:ref:`outfiles`):

.. code-block:: c

   sam_header_inputs_t inputs;
   sam_header_inputs_init(&inputs);
   sam_header_inputs_add_files(&inputs, files, n_files, logger);
   char *cl = sam_header_command_line(argc, argv);
   for (size_t d = 0; d < n_destinations; d++) {
     strbuf_clear(&text);
     strbuf_append_str(&text, "@HD\tVN:1.6\tSO:unknown\n");
     sam_header_append_read_groups(&text, &inputs, names[d]);
     sam_header_append_programs(&text, &inputs, "demux", VERSION, cl);
     strbuf_append_char(&text, '\0');        /* 0-terminated */
     strbuf_clear(&bam);
     ubam_writer_header(&bam, text.s);
     strbuf_clear(&compressed);
     bgzf_compressor_run(bgzf, bam.s, bam.len, &compressed);
     outfiles_set_bam(out, d, compressed.s, compressed.len);
   }
   free(cl);
   sam_header_inputs_free(&inputs);

With an input ``in.bam`` with the read groups ``run1_A`` and ``run1_B``
and the basecaller as last program, the header of the sample ``S1``
written by ``demux -o out in.bam`` is (fields separated by tabs)::

   @HD  VN:1.6  SO:unknown
   @RG  ID:S1  SM:S1
   @RG  ID:S1.run1_A  PL:ONT  SM:S1
   @RG  ID:S1.run1_B  PL:ONT  SM:S1
   @PG  ID:basecaller  PN:dorado  VN:0.9
   @PG  ID:demux  PN:demux  PP:basecaller  VN:1.0  CL:demux -o out in.bam

SAM tags in read titles (``sam-tags.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The titles of the records (FASTQ titles, the titles returned by the BAM
reader, and those given to the uBAM writer) contain the read name and SAM
tags in text format, separated by tabs (as in the FASTQ of Dorado, for
``minimap2 -y`` or ``samtools import -T '*'``). Functions for building
and checking them:

``bool sam_tag_is_valid(const char *field, size_t len)``
   whether a field is a tag in SAM text format: a two-character name (a
   letter and a letter or digit), a type (``A``, ``i``, ``f``, ``Z``,
   ``H`` or ``B``) and the value, separated by ``:`` (e.g.
   ``RG:Z:run1``);
``void sam_tag_append_string(strbuf_t *out, const char *tag, const char *value, size_t len)``
   append a tab and the string tag ``<tag>:Z:<value>``;
``void sam_tag_append_int(strbuf_t *out, const char *tag, long value)``
   append a tab and the integer tag ``<tag>:i:<value>``;
``void sam_tag_start(strbuf_t *out, const char *tag, char type)``
   append a tab and the start of a tag (e.g. ``\tQT:Z:``), for values
   built piece by piece (e.g. the qualities of several indices).

For example, a demultiplexer can append to the title of a read its
sample barcode and the number of a molecule:

.. code-block:: c

   strbuf_append(&title, name, name_len);
   sam_tag_append_string(&title, "BC", index, index_len);
   sam_tag_append_int(&title, "mk", molecule);

