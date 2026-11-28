Troubleshooting
===============

Each message of flucut, its cause and what to do. Errors end the program
with exit status 1.

Errors in the command line
--------------------------

``unknown scheme '<name>'``
   The name given to ``-s`` is not a built-in scheme; the help lists
   them.

``'primer' or 'terminus' expected``
   The value of ``-T`` is neither.

``a fraction below 1, a whole number of bases or auto expected``
   The value of ``-e`` is negative, not a number nor ``auto``, or a
   fraction above 1 (e.g. 1.5).

``argument must be an integer ...``
   A value of ``-w``, ``-m`` or ``-j`` out of its range.

``Option '--probe-len' / '--probe-errors' / '--terminal-window' ...``
   A value out of its range (see :doc:`cli_reference`).

``the output directory is empty``
   ``-o ""`` was given.

``Option '--no-platform-check': a platform is required``
   ``--no-platform-check`` turns off the check of a platform given; with
   ``-p auto`` (the default) there is none to check. Add ``-p nanopore``
   or ``-p illumina``.

Errors in the inputs
--------------------

``<file>: not a file of reads``
   An input which is not a directory must end in ``.fastq``, ``.fq``
   (optionally ``.gz``) or ``.bam``.

``No FASTQ or BAM files in the input group <name>`` (warning)
   A directory given as input, or one of its barcode directories,
   contains no file with these suffixes: its output is empty.

``The input directory <dir> contains barcode directories and <n> FASTQ or BAM file(s) outside them`` (warning)
   Only the barcode directories are trimmed; the files beside them are
   not read.

``The input group name <name> is not unique``
   Two groups have the same name (e.g. two directories ``barcode01`` in
   different runs, or a file and a directory with the same name): trim
   them in separate runs, with different output directories.

``<file>: the input would be overwritten by the output``
   The output directory is the directory of an input with the same
   name: choose another output directory.

Errors in the references file
   ``the sequence of "<name>" is empty``, ``no sequence found``,
   ``invalid base``, ``sequence before the first header line``: the
   file given to ``-r`` is not a valid FASTA file.

``could not download the built-in references ...`` (warning)
   The default reference set could not be downloaded (no network, or
   ``curl`` is missing) and is not cached yet, so the run proceeds without
   references (every primer-derived base removed). Connect to the network
   for the first run (it is cached afterwards, under
   ``$XDG_CACHE_HOME/flucut`` or ``~/.cache/flucut``), pass your own with
   ``-r <file>``, or ``-r off`` to turn references off and silence this.

``the references look like coding sequences without the UTRs ...`` (warning)
   The references given to ``-r`` start at a start codon and end at a
   stop codon: they are coding sequences (CDS), without the untranslated
   regions (UTRs) of the segments. flucut locates a read by the sequence
   just after its 5' primer, which in a complete segment is the 5' UTR;
   a CDS reference does not contain it, so the full-length (terminal)
   reads cannot be located (they are counted as not located) and only
   the internally primed reads, whose probe falls inside the ORF, map —
   the internal-priming classification and the coverage plot are then
   misleading. Use **complete segment references, with the UTRs**:
   ``scripts/fetch_references.py`` downloads them, and the records of a
   complete segment begin with the conserved terminus (``AGCRAAAGCAGG``),
   not with ``ATG``.

Messages about the reads
------------------------

``the reads look like Illumina reads, not nanopore (-p nanopore): ...``
   ``-p nanopore`` was given, but the first reads are of a fixed length
   no longer than 320 bases, as an Illumina run gives
   (:ref:`platform-check`). flucut stops before trimming. Leave out
   ``-p`` or pass ``-p illumina``.

``the reads look like nanopore reads, not Illumina (-p illumina): ...``
   ``-p illumina`` was given, but at least 1% of the first reads are
   longer than 600 bases, which no Illumina read is. flucut stops before
   trimming. Leave out ``-p`` or pass ``-p nanopore``.

``the reads may not be <platform> reads: ...`` (warning)
   The first reads lean to the other platform without being clearly of
   it: e.g. no read longer than 600 bases but no fixed length (Illumina
   reads trimmed to many lengths, or joined pairs; nanopore reads of a
   library of short molecules), or a few long reads among Illumina-like
   ones. flucut runs with the platform given (with ``-p auto``,
   nanopore); the reports carry the remark. If the reads are of the
   other platform, run again with its ``-p``.

``Platform: Illumina, read from the reads (-p auto): ...``
   Not a problem: the reads were read as Illumina reads, and trimmed
   with the Illumina settings (window 10, minimum length 1, no adapter
   search, 15% errors in a primer with ``-e auto``).

``not located in the references: <n> left, <n> right``
   Read ends whose probe was not found in the references. A few are
   expected (reads with many errors, dimers); many mean that the
   references do not fit the run (e.g. a missing subtype): download
   references of the subtypes of the run.

``<reference>: <n> of <m> complete amplicons shorter than 90% of the reference``
   Complete amplicons much shorter than their segment: internal priming
   or deletions (:doc:`internal_priming`). Look at the reads of
   ``lengths.tsv`` and at the ``ml``/``mr`` tags.

Many reads without primers (``no_end`` in ``summary.tsv``)
   The scheme does not match the primers of the run (``-s``), the reads
   were already trimmed (e.g. by the basecaller), or the error rate of
   the reads is higher than ``-e`` allows (with ``-e auto``, see the rate
   measured in the report).

``only N% of the reads placed on a segment`` (remark on a sample, tinted red; column *placed* of the overview)
   Few of the reads could be placed on a segment of the references, so
   most of them are of something else: another organism, nucleic acid
   of the host, or an amplification which did not give the expected
   amplicons. What the report says about the segments, the coverage and
   the subtype of such a sample rests on the few reads that were placed.
   The reads with a primer at both ends which were not placed are
   counted as ``unplaced`` (:ref:`reads-by-category`) and are written to
   ``<name>.unplaced.fastq.gz`` by ``--by-category``, so
   ``<name>.full.fastq.gz`` holds only reads placed on a segment. If the sample
   should be influenza, check that the references cover its subtype and
   its type (``--virus``).
