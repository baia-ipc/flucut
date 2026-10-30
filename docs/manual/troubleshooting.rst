Troubleshooting
===============

Each message of flucut, its cause and what to do. Errors end the program
with exit status 1.

Errors in the command line
--------------------------

``--trim terminus requires --references``
   Keeping the termini without references would keep the primer bases
   of the internally primed reads (the artefact). Give the references
   (``-r``), or use the default ``--trim primer``.

``unknown scheme '<name>'``
   The name given to ``-s`` is not a built-in scheme; the help lists
   them.

``'primer' or 'terminus' expected``
   The value of ``-T`` is neither.

``a fraction below 1 or a whole number of bases expected``
   The value of ``-e`` is negative, not a number, or a fraction above 1
   (e.g. 1.5).

``argument must be an integer ...``
   A value of ``-w``, ``-m`` or ``-j`` out of its range.

``Option '--probe-len' / '--probe-errors' / '--terminal-window' ...``
   A value out of its range (see :doc:`cli_reference`).

``the output directory is empty``
   ``-o ""`` was given.

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

Messages about the reads
------------------------

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
   the reads is higher than ``-e`` allows.
