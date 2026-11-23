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
(:doc:`reports`): ``flucut_report.html``, the HTML report,
``flucut_dashboard.html``, the dashboard with which the run can be
followed, and ``flucut_summary.json``, the same content for programs,
with the statistics of all reads and of each group.

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
``a5``   the kit of the adapter removed at the start of the read
         (``nb``, ``rbk`` or ``lsk``; only at an end without a
         primer, see :doc:`how_it_works`)
``a3``   the same, at the end of the read
``ol``   with ``-r``: origin of the left end (``terminal``,
         ``internal`` or ``unknown``)
``or``   the same for the right end
``ml``   with ``-r``: where the probe of the left end mapped,
         ``<reference>:<position>:<strand>`` (1-based position on the
         plus strand)
``mr``   the same for the right end
``ct``   with ``-r``: the category of the read: ``full``,
         ``internal``, ``deletion`` or ``other`` (see
         :ref:`reads-by-category`)
``sp``   ``<k>/<n>``: the record is part ``k`` of the ``n`` parts of a
         read which contained several amplicons (see below)
======== ============================================================

The left and right ends are those of the amplicon: for a read of the
minus strand, the left end is at the end of the read. IRMA accepts the
tags (:doc:`irma`).

A read which contains several amplicons, joined end to end, is split at
their junctions (:doc:`internal_priming`), and each part is written as a
record of its own. The name of the read (the first word of the title) is
followed by ``_1``, ``_2``, ... so that the names stay unique, and the
tags, ``tl`` and ``tr`` included, describe the part. ``--no-split`` writes
one record for each read.

.. _reads-by-category:

Reads by category
-----------------

With references, each read written belongs to one of four categories,
given by its ``ct`` tag. They are the groups of the coverage plot of the
HTML report (:doc:`reports`), which describes how the reads are told
apart:

============ ========================================== ================
``ct``       reads                                      in the plot
============ ========================================== ================
``full``     full templates: the amplicon started at a  terminal +/-
             segment end                                (green)
``internal`` internal priming: the amplicon started     internal +/-
             inside the segment                         (amber, blue)
``deletion`` internal deletions: the two ends of a      internal
             segment joined, e.g. defective             deletion
             interfering genomes                        (purple)
``other``    reads without a 5' primer (fragments of    grey, or not
             unknown origin), and reads which could     shown
             not be located on the references
============ ========================================== ================

The category describes where the amplicon started (its 5' end), on
either strand; the strand is in the ``fs`` tag.

``--by-category`` writes the reads of each input to four files, one for
each category, instead of ``<name>.fastq.gz``::

   flucut --by-category -o trimmed barcode01/

   trimmed/barcode01.full.fastq.gz
   trimmed/barcode01.internal.fastq.gz
   trimmed/barcode01.deletion.fastq.gz
   trimmed/barcode01.other.fastq.gz

All four files are always written; one without reads is an empty gzip
file. In each file the reads are in the order of the input, and together
the four files contain the same reads as ``<name>.fastq.gz`` would. The
parts of a split read are each written to the file of their own
category. The tables and the reports have one part for each input, not
for each file; the reads of each file are counted in the columns
``out_full``, ``out_internal``, ``out_deletion`` and ``out_other`` of
``summary.tsv``.

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
   * - ``split``
     - reads split into several amplicons; ``reads_in`` and ``bases_in``
       count the reads of the input, all the other columns count each
       part as a read
   * - ``left_terminal``, ``left_internal``, ``left_unknown``, ``right_terminal``, ``right_internal``, ``right_unknown``
     - with ``-r``: the left and right ends found, by origin; ``NA`` without references
   * - ``short_complete``
     - with ``-r``: complete amplicons shorter than 90% of their reference; ``NA`` without references
   * - ``out_full``, ``out_internal``, ``out_deletion``, ``out_other``
     - with ``-r``: the reads written, by category (:ref:`reads-by-category`); their sum is ``reads_out``; ``NA`` without references

lengths.tsv
-----------

The lengths of the complete amplicons (both primers found; the template
between the two primers), in bins of 10 bases, tab-separated, with a
header: ``output``, ``group`` (the reference to which the probe of the
left end, or else of the right end, mapped; ``not_located``; or ``all``
without references), ``length`` (the start of the bin) and ``reads``.
Only the bins with reads are written; the last bin collects the lengths
from 4000 bases on.

flucut_report.html, flucut_dashboard.html, flucut_summary.json
---------------------------------------------------------------

The reports of the run, always written (:doc:`reports`): the report at
the end of the run; the dashboard and the JSON summary also while it is
in progress, if it has several groups. They contain the settings of the
run and, for each group (a sample), its state and subtype, the reads by
primers found and by origin, their lengths and qualities before and
after the trimming, the segments with their depth, coverage by read
origin and complete amplicons, and the internal priming sites. Each is
written to a temporary file (suffix ``.tmp``) and then renamed.

Messages
--------

The messages go to the standard error (or to the file of ``-L``): the
counts of each input, the internal priming found, and the length check
(:doc:`internal_priming`). In a terminal, the banner of flucut is shown
first (not with ``-q`` or ``-L``).
