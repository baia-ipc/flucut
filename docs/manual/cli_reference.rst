Command line reference
======================

::

   flucut [options] <input>...

Arguments
---------

``<input>``
   A FASTQ file (optionally gzipped), an unaligned BAM file, or a
   directory of such files (recognised by the suffixes ``.fastq``,
   ``.fq``, optionally followed by ``.gz``, and ``.bam``; e.g. a barcode
   directory of a MinKNOW run, whose files are trimmed into one output).
   A directory with subdirectories whose name starts with ``barcode``
   (e.g. the ``fastq_pass`` directory of a MinKNOW run) gives one output
   for each of them, and the files beside them are not read. Each group
   of files gives one output (:doc:`output`).

Options
-------

``-o, --outdir <dir>``
   The output directory, created if needed. Default: ``flucut-output``.

``-s, --scheme <name>``
   The primer scheme (:doc:`primers`): ``uni12-mbt13`` (default),
   ``mbt`` or ``uni``.

``-r, --references <file>``
   A FASTA file of segment references: each read end is classified as
   terminal or internally primed (:doc:`internal_priming`), and the
   trimming mode follows from this. **With** ``-r`` the conserved
   termini are kept at the ends confidently classified as terminal
   (complete segment ends) and the whole primer is removed elsewhere
   (an internal or an unclassified end is cut when unsure). **Without**
   ``-r`` every primer-derived base is removed. There is no separate
   option to choose the mode (:doc:`how_it_works`).
   The messages and reports name a reference with its segment when its
   description (the header after the name) says it: ``segment <n>``
   (as in the NCBI records), or a gene or protein name (e.g. ``HA``,
   ``hemagglutinin``, or the IRMA names like ``A_HA_H5``). Default:
   none.

``--adapters <kits>``
   The sequencing adapters of Oxford Nanopore searched at the read ends
   at which no primer is found (:doc:`how_it_works`): a comma-separated
   list of ``nb`` (native barcoding), ``rbk`` (rapid barcoding) and
   ``lsk`` (ligation, without barcodes), or ``none``. Default:
   ``nb,rbk,lsk``.

``-e, --max-errors <x>``
   The errors allowed in a primer: a whole number of bases, below 1 a
   fraction of the length of the primer, or ``auto``: decided from the
   error rate of the primer alignments in the first 1000 reads of the
   first input file (with 0.2): 0.2 (2 errors in Uni12), or 0.3 (3
   errors) if that rate is at least 4% (:doc:`how_it_works`). The value
   decided and the rate measured are in the messages (``-v``) and in the
   reports. Default: ``auto``.

``-p, --platform <p>``
   ``nanopore`` (default): reads of whole amplicons, with flanking bases
   and adapters before the primers. ``illumina``: reads of fragments of
   the amplicons (e.g. after tagmentation), whose adapters were removed
   by the instrument: a primer, when there is one, is at the start of
   the read; the window is 10 bases and the Nanopore adapters are not
   searched (unless ``-w`` or ``--adapters`` are given); with
   ``--references``, the read ends are located by their k-mers only,
   without the exhaustive search of those the k-mers do not locate (with
   the few errors of Illumina reads, a read end of a segment shares
   exact k-mers with its reference; the reads of the host, which belong
   to no reference, would each be aligned to every reference: 3 times
   slower on a run of milk samples). With the
   window of 100 bases, a primer-like site in the first 100 bases of an
   Illumina read would be taken for the primer (:doc:`performance`).

``-w, --window <n>``
   A primer is accepted at an end only if it begins (left end) or ends
   (right end) within ``<n>`` bases of that end of the read. Default:
   100; with ``-p illumina``, 10.

``-m, --min-length <n>``
   Reads shorter than ``<n>`` bases after the trimming are not written.
   Default: 1 (empty reads, e.g. primer dimers, are never written).

``-j, --threads <n>``
   The number of worker threads (1 to 256; default 4), each trimming
   batches of 1000 reads (:doc:`how_it_works`). The output does not
   depend on it.
   Default: 4.

``--compress-level <n>``
   The gzip compression level of the outputs, 1 to 12 (the levels of
   libdeflate, with which they are compressed). A higher level gives
   smaller files, and takes longer (:doc:`performance`). Default: 6.

``--no-tags``
   Do not append the tags to the read titles.

Classification against the references
-------------------------------------

These options are used with ``--references``.

``--probe-len <n>``
   The number of bases after the primer of an end located in the
   references. At least 16. Default: 32.

``--probe-errors <x>``
   The fraction of the probe allowed to be errors (0 to below 1).
   Default: 0.2.

``--terminal-window <n>``
   A probe mapped within ``<n>`` bases of an end of a segment is
   terminal. Default: 60.

``--exhaustive-search``
   Align each probe to the whole of every reference, without the k-mer
   index (slower; for checking the index).

Live mode
---------

The live mode (:doc:`live`): its options can only be given with
``--live``.

``--live``
   Trim the files of the input directories as they are written (e.g. by
   MinKNOW), appending to the outputs and updating the tables and
   reports after each file; stop at the end of the run (a file
   ``final_summary_*.txt`` in an input directory or in its parent),
   after ``--timeout``, or with Ctrl-C. The inputs must be directories.

``--timeout <t>``
   Also stop after ``<t>`` without new files (e.g. ``30s``, ``60m``,
   ``2h``; a number: minutes). Default: none.

``--min-age <t>``
   Trim a file only ``<t>`` after its last change (e.g. ``30s``,
   ``2m``; a number: seconds). Default: 60 seconds.

``--scan-interval <s>``
   Seconds between the scans of the input directories (at most 3600).
   Default: 10.

``--resume``
   Skip the files already trimmed by an interrupted live run with the
   same output directory, continuing its outputs and counts.

Other options
-------------

``-q, --quiet``
   Only error messages.
``-v, --verbose``
   Additional information (e.g. the values of the options).
``-D, --debug``
   Debugging information.
``-L, --logfile <file>``
   The messages to a file instead of the standard error.
``-h, --help``
   The help, with the list of the schemes.
``-V, --version``
   The version.

Only one of ``-q``, ``-v`` and ``-D`` can be given. Invalid values of the
options are errors (exit status 1), reported with the usage.

The messages go to the standard error; in a terminal they are colored by
their kind:

- ``SUCCESS`` (green): the reads written for each output, and the
  reports at the end;
- ``WARNING`` (magenta): something to look at, e.g. complete amplicons
  shorter than the closest reference (below 90% of its length: internal
  priming or deletions; with ``-r``), named by their segment when the
  reference says it, as in ``HA segment: 4 of 9 complete amplicons
  shorter than the closest reference (NC_004908.1, 1714 bases)``, or
  files of reads ignored beside the barcode directories;
- ``ERROR`` (red): the program stops;
- with ``-v``, also ``PARAM`` (blue: the values of the options),
  ``PHASE`` (bright yellow: the steps, e.g. the trimming of each file)
  and ``INFO`` (yellow: the scheme and its primers, the input groups,
  the references read, the internally primed ends, the lengths of the
  amplicons without references, the tables written).

Without ``-v``, a run thus prints one green line for each output, a
warning if any, and a green line for the reports.
