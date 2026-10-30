Outputs
=======

The inputs are read as groups of files, each of which gives one
output, ``<outdir>/<name>.fastq.gz``:

* a file of reads given directly is a group, named as the file without
  its suffix (``.fastq``, ``.fq``, optionally ``.gz``, or ``.bam``);
* a directory is a group of its files of reads, named as the directory;
* unless the directory contains subdirectories whose name starts with
  ``barcode`` (the ``fastq_pass`` directory of a MinKNOW run, with its
  ``barcodeNN`` directories): then each of them is a group, named as the
  subdirectory; the files beside them are not read (a warning is given),
  and the other subdirectories (e.g. ``unclassified``) are ignored.

The files of a directory are read in the order of their names. Two
groups with the same name are an error, as one would overwrite the
other. ``summary.tsv`` and ``lengths.tsv`` have one part for each group.

The output directory also contains the reports of the run
(:doc:`reports`): ``flucut_report.html``, the HTML report, and
``flucut_summary.json``, the same content for programs, with the
statistics of all reads and of each group.

The trimmed reads
-----------------

``<name>.fastq.gz`` (gzip, FASTQ) contains the reads which are not empty
and not shorter than ``-m`` bases after the trimming, in the order of the
input, in their original orientation, with their quality values cut as
the sequence. Their title is the original title followed by tags in the
SAM text format, separated by tabs (as in the FASTQ written by Dorado),
unless ``--no-tags`` is given:

======== ============================================================
tag      meaning
======== ============================================================
``fs``   strand of the amplicon: ``+``, ``-``, or ``.`` (no primer)
``fl``   the primer found at the left end (Uni12 side)
``fr``   the primer found at the right end (Uni13 side)
``tl``   bases removed at the start of the read
``tr``   bases removed at the end of the read
``ni``   primer matches inside the read (if any)
``ol``   with ``-r``: origin of the left end (``terminal``,
         ``internal`` or ``unknown``)
``or``   the same for the right end
``ml``   with ``-r``: where the probe of the left end mapped,
         ``<reference>:<position>:<strand>`` (1-based position on the
         plus strand)
``mr``   the same for the right end
======== ============================================================

The left and right ends are those of the amplicon: for a read of the
minus strand, the left end is at the end of the read. IRMA accepts the
tags (:doc:`irma`).

summary.tsv
-----------

One line for each input, tab-separated, with a header:

.. list-table::
   :header-rows: 1
   :widths: 38 62

   * - column
     - meaning
   * - ``output``
     - the name of the input
   * - ``reads_in``
     - reads read
   * - ``bases_in``
     - their bases
   * - ``reads_out``
     - reads written
   * - ``bases_out``
     - their bases
   * - ``both_ends``
     - reads with the primers of both ends
   * - ``one_end``
     - reads with the primer of one end
   * - ``no_end``
     - reads without primers
   * - ``minus_strand``
     - reads found on the minus strand of the amplicon
   * - ``dimers``
     - primer dimers (both primers, nothing between)
   * - ``too_short``
     - reads not written: empty or shorter than ``-m`` (dimers included)
   * - ``with_internal_match``
     - reads with primer matches inside them
   * - ``left_terminal``, ``left_internal``, ``left_unknown``, ``right_terminal``, ``right_internal``, ``right_unknown``
     - with ``-r``: the left and right ends found, by origin; ``NA`` without references
   * - ``short_complete``
     - with ``-r``: complete amplicons shorter than 90% of their reference; ``NA`` without references

lengths.tsv
-----------

The lengths of the complete amplicons (both primers found; the template
between the two primers), in bins of 10 bases, tab-separated, with a
header: ``output``, ``group`` (the reference to which the probe of the
left end, or else of the right end, mapped; ``not_located``; or ``all``
without references), ``length`` (the start of the bin) and ``reads``.
Only the bins with reads are written; the last bin collects the lengths
from 4000 bases on.

flucut_report.html, flucut_summary.json
---------------------------------------

The reports of the run, always written, at its end: the run
parameters and, for all reads and for each group, the reads at each
step, the primers found at each end, the internal priming sites (with
references) and the complete amplicons, with the histograms of their
lengths (:doc:`reports`). Each is written to a temporary file (suffix
``.tmp``) and then renamed.

Messages
--------

The messages go to the standard error (or to the file of ``-L``): the
counts of each input, the internal priming found, and the length check
(:doc:`internal_priming`). In a terminal, the banner of flucut is shown
first (not with ``-q`` or ``-L``).
