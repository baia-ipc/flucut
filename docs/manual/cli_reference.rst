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
   Each input gives one output (:doc:`output`).

Options
-------

``-o, --outdir <dir>``
   The output directory, created if needed. Default: ``flucut-output``.

``-s, --scheme <name>``
   The primer scheme (:doc:`primers`): ``uni12-mbt13`` (default),
   ``mbt`` or ``uni``.

``-T, --trim <policy>``
   ``primer`` (default): remove every primer-derived base. ``terminus``:
   keep the conserved termini; requires ``--references``
   (:doc:`how_it_works`).

``-r, --references <file>``
   A FASTA file of segment references: the read ends are classified as
   terminal or internally primed (:doc:`internal_priming`). With
   ``--trim primer`` the classification is only reported; with
   ``--trim terminus`` it decides the trimming of the internal ends.
   Default: none.

``-e, --max-errors <x>``
   The errors allowed in a primer: a whole number of bases, or, below 1,
   a fraction of the length of the primer. Default: 0.2.

``-w, --window <n>``
   A primer is accepted at an end only if it begins (left end) or ends
   (right end) within ``<n>`` bases of that end of the read. Default:
   100.

``-m, --min-length <n>``
   Reads shorter than ``<n>`` bases after the trimming are not written.
   Default: 1 (empty reads, e.g. primer dimers, are never written).

``-j, --threads <n>``
   The number of threads (1 to 256). The output does not depend on it.
   Default: 4.

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
