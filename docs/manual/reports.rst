.. _reports:

Reports
=======

Besides the trimmed reads, the tables and the messages, flucut always
writes three files to the output directory:

* ``flucut_dashboard.html``, a dashboard to follow the run;
* ``flucut_report.html``, the report of the run, written at its end;
* ``flucut_summary.json``, the same information for programs
  (pipelines, other dashboards).

Each is written to a temporary file which is then renamed, so that a
partial file is never read. Both pages are single, self-contained files:
the plots are inline SVG and no external resources are used, so that
they can be opened offline in any browser, archived, and sent by
e-mail. They follow the light or dark theme of the browser.

The dashboard
-------------

A run with several samples (several input files or directories, or the
barcode directories of a run) can be followed while it is in progress.
The dashboard is written before the first sample and after each one (in
the :doc:`live mode <live>`, after each input file), and flucut says
where it is as soon as it exists::

   SUCCESS | Dashboard: out/flucut_dashboard.html (open it in a web
     browser to follow the run)

It reloads itself in the browser and shows one thing at a time: the
**overview** of the run, the **HA** in all samples, or one **sample**,
chosen with the links at the top or by the name of a sample in the
table of the overview. The choice is kept in the page address, so it
survives the reloads. Samples still to be processed are shown as
waiting. At the end of the run the dashboard stays, in its last state,
with a link to the report (the JSON summary has ``"running": true``
until then).

.. figure:: figures/dashboard.*
   :width: 100%

   The dashboard of a run of four barcodes (simulated reads), at the
   end of the run.

The report
----------

The report has the same parts as the dashboard, all on one page: the
run overview, the HA in all samples, then each sample. With one sample
only, it is that sample.

The figures of this chapter show a simulated run of four barcodes
(``docs/manual/figures/README.md`` says how it was made): two complete
samples, one with few reads of a polymerase segment, and one with too
few reads.

How to read it
~~~~~~~~~~~~~~

A few conventions hold everywhere:

* **Segments are judged by their depth.** A segment is *complete* when
  at least 90% of its positions have a depth of 50 or more (what
  FLU-minion of IRMA needs to support a base), *partial* from 50%, and
  *low depth* below. Depth and coverage are estimated from where flucut
  locates each read on its reference (its start and its length), not
  from a base-by-base mapping.
* **A sample is summarised by a count**: its complete segments, of the
  eight (e.g. ``7/8``), green when all are complete, amber when some
  are partial or low, red when some are missing. The count is grey when
  the segments cannot be judged: fewer than 1,000 reads, or half or
  more of the read ends not matched to the references (e.g. a sample
  which is not influenza A). In the live mode it is grey as long as
  reads are still being added.
* **Reads have one of four origins**, always in the same colours: *full
  template* (blue: the amplicon started at an end of the segment),
  *internally primed* (magenta: it started inside the segment), *internal
  deletion* (violet: the two ends of the read lie much farther apart on
  the segment than the read is long, as in defective interfering
  genomes) and *other* (grey: no 5' primer, or not located). They are
  the categories of the ``ct`` tag and of ``--by-category``
  (:ref:`reads-by-category`).
* **Every percentage names what it is a share of** (of the reads after
  splitting, of the depth of a segment, of the read ends).
* ``0.0%`` and ``100.0%`` mean exactly none and exactly all: a share
  which is not zero, but below 0.1%, is written ``<0.1%``, and one above
  99.9%, but not all, ``>99.9%``.

Details which are rarely needed (all the settings, the smaller internal
priming sites, the read ends by origin, the terms) are folded: click
their line to open them.

The head
~~~~~~~~

The platform, the primer scheme, the trimming policy, the number of
references and the version on one line; under it, when the titles of
the reads give them, the name of the experiment, the flow cell and the
basecalling model. *Settings* opens the full table: command line,
references, errors allowed in a primer (and how ``-e auto`` decided),
window, minimum length, probe, threads.

Run overview
~~~~~~~~~~~~

.. figure:: figures/report_overview.*
   :width: 100%

   The overview of the simulated run.

**The samples by state**, on one line (complete, partial, with segments
missing, few reads, not matched to the references, waiting).

**Tiles** with the throughput of the run: samples, subtype (and how
many samples have each, if several), reads and bases, median read
length (with the N50 and the median quality), reads written, reads with
both primers.

**The table of the samples**, one row each; the name links to the
sample:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - column
     - content
   * - segments
     - the complete segments, of those expected
   * - subtype
     - the H of the reference closest to its HA and the N of that
       closest to its NA; with a question mark when the HA or the NA
       has a low depth, or the sample cannot be judged
   * - reads, N50, quality
     - the reads of the input (the bar compares the samples), the N50
       of their lengths and the median of their qualities
   * - PB2 ... NS
     - for each segment, the percentage of its positions at a depth of
       50 or more, tinted by state (green: complete, amber: partial,
       red: low depth or not found; no tint for a sample which cannot
       be judged)
   * - int. primed, int. deletions
     - the segment with the greatest share of its depth from internally
       primed reads, from reads with an internal deletion; bold from
       10% and from 50%. Pointing at the cell lists all the segments.
       Not given for a segment with a mean depth below 20.
   * - reads by origin
     - the reads written, by origin, as a small bar
   * - both primers
     - the reads with the primers of both ends

**Median depth by segment**: a dot for each sample over each segment,
on a logarithmic axis, with a line at 50. A segment which is low in all
samples points to the amplification; one low in a single sample points
to that sample. The dots of the samples which cannot be judged are
hollow.

Nothing else is pooled across samples: coverage, internal priming sites
and amplicon lengths are those of each sample.

HA in all samples
~~~~~~~~~~~~~~~~~

.. figure:: figures/report_ha.*
   :width: 100%

   The HA of the first two samples of the simulated run: 36% and 8% of
   the depth come from reads internally primed at position 1090.

The HA decides the subtype and is where internal priming matters most
(:doc:`internal_priming`), so it has a view of its own: for each sample,
its subtype, the strip of its HA (see *Segments* below), side by side
for comparison. The view is left out when the run or the report is
limited to another segment (``--segment``, ``--report-segment``).

A sample
~~~~~~~~

.. figure:: figures/report_sample.*
   :width: 100%

   The top of a sample: state and remarks, tiles, reads, internal
   priming sites, bases by segment.

**State and remarks.** The count of complete segments and, when they
apply, short remarks:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - remark
     - shown when
   * - few reads
     - fewer than 1,000 reads
   * - no primer in N% of the reads
     - 20% or more: a scheme which does not match the primers of the
       run (``-s``), reads already trimmed, or more errors than ``-e``
       allows
   * - N% of the read ends not located in the references
     - 50% or more
   * - HA 36% of depth internally primed
     - 10% or more of the depth of a segment
   * - PB2 90% of depth from reads with an internal deletion
     - 50% or more of the depth of a segment

**Tiles**: the subtype with the references closest to the HA and to the
NA, the reads and bases, the median read length, the reads written, the
reads with both primers.

**Reads.** Two bars: the reads by primers found (both, one, none) and
by origin (with the reads not written: too short after the trimming,
dimers included). Under them the smaller counts: the primer of each end
by name, minus strand, reads split into several amplicons, primer
dimers, too short, primer matches inside the read, adapters removed.
*Read ends by origin* opens the same division for the ends with a
primer (terminal, internal, not located).

**Internal priming sites.** The sites with 10 read ends or more, the
largest first (at most five; the others are folded): segment, end
(left: Uni12, right: Uni13), position of the first base after the
primer (1-based on the plus strand of the reference, as in the
``ml``/``mr`` tags; read ends at most 5 bases apart are one site), reads
and their share of all reads.

**Bases by segment.** The share of the located bases in each segment:
the balance of the amplification.

.. figure:: figures/report_segments.*
   :width: 100%

   The first segments of a sample.

**Segments.** One strip for each segment (a segment without located
reads is *not found*):

* at the left, its state, its reference, and its numbers: median depth,
  share of the positions at a depth of 50 or more, share covered, and
  the shares of its depth from internally primed reads, from reads with
  an internal deletion and from reads of the minus strand;
* in the middle, the **coverage by read origin**: the depth along the
  segment, stacked by origin. The reads of unknown origin and the full
  templates are at the bottom, as flat bands; the internally primed
  reads are above them, so that internal priming shows as a step up at
  its site, which is marked with its position when 5% or more of the
  depth is internally primed; reads with an internal deletion cover the
  two ends only. The dashed line is the depth of 50, and the stretches
  below it are marked in red under the axis;
* at the right, the **complete amplicons by length**, with the length of
  the segment marked: amplicons shorter than 90% of it (grey) come from
  internal priming or from deletions.

The coverage needs references which are **complete segments, with the
UTRs**: a read is located by the sequence just after its 5' primer,
which in a full segment is the 5' UTR. Against coding sequences only
the internally primed reads map, and flucut warns
(:doc:`troubleshooting`). How the reads are placed is described in
:doc:`how_it_works`.

.. figure:: figures/report_reads.*
   :width: 100%

   The lengths and qualities of the reads of a sample.

**Read length and quality.** A table of the reads of the input and of
the reads written: reads, bases, median and mean length, N50, shortest
and longest, median quality, and the reads with a quality of 10 and of
20 or more; and the distributions of the lengths and of the qualities,
the input as bars and the reads written as an outline. The quality of a
read is the Phred score of the mean error probability of its bases, as
in MinKNOW. Reads are shorter after the trimming, and usually of a
higher quality, as adapters and primers are removed. Under the table,
what the titles of the reads say about the sequencing: flow cell, run,
basecalling model (which names the pore, the chemistry, the speed and
the accuracy, e.g. ``dna_r10.4.1_e8.2_400bps_sup@v4.3.0``), sample name
and experiment as given in MinKNOW, and the time of the first read; for
Illumina reads, the instrument, the run number and the flow cell. A
value which differs among the reads is given as *several*.

Without references
~~~~~~~~~~~~~~~~~~

With ``-r off`` the read ends are not classified: there is no subtype,
no segment and no internal priming site. The report has the reads by
primers found, the reads written, the lengths of all complete amplicons
(internal priming shows there as amplicons much shorter than their
segment) and the read lengths and qualities.

One segment
~~~~~~~~~~~

With ``--report-segment`` (e.g. ``--report-segment HA``) the parts of
the report which list the references show one segment only; the counts
of the reads are unchanged. With ``--segment``, which also limits the
classification to the references of one segment, only that segment is
expected: a sample is then ``1/1`` when it is complete
(:doc:`cli_reference`).

The JSON summary
----------------

``flucut_summary.json`` contains the same information as the HTML
report, in a documented format (version 1, described in
``docs/json_summary.md`` of the repository). Its main keys are:

* ``format``, ``program``, ``version``, ``updated``: the format
  version, the program and the time of writing;
* ``command``, ``scheme``, ``trim_policy``, ``references``,
  ``settings``: the parameters of the run (among them ``segment`` and
  ``report_segment``, the segment of ``--segment`` and of
  ``--report-segment``, ``null`` for all segments);
* ``all``: the statistics of all reads, and ``groups``: those of each
  group, with ``reads`` (the reads at each step), ``read_lengths`` and
  ``read_quality`` (of the input and of the reads written),
  ``sequencing`` (the metadata of the run in the titles of the reads),
  ``sample`` (with references: state, subtype, remarks and the eight
  segments), ``primers``, ``ends`` and ``internal_sites`` (with
  references), ``short_complete`` and ``complete_amplicons`` (with the
  histograms of the lengths).

For example, the internal priming sites of each group with ``jq``::

   jq -r '.groups[] | .name as $g | .internal_sites[]? |
          [$g, .reference, .end, .position, .reads] | @tsv' \
       out/flucut_summary.json

which gives for the simulated reads::

   reads  KF918483.1   left  1090  224
   reads  KF918483.1   left  1158  6
   reads  NC_002023.1  left  1595  1
