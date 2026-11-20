.. _reports:

Reports
=======

Besides the trimmed reads, the tables and the messages, flucut always
writes two reports to the output directory:

* ``flucut_report.html``, a report for people;
* ``flucut_summary.json``, the same information for programs
  (pipelines, dashboards).

Both are written at the end of the run, to a temporary file which is
then renamed, so that a partial report is never read. A run with several
samples (several input files or directories, or the barcode directories
of a run) can be followed while it is in progress: the tables and the
reports are also written before the first sample and after each one.
Until the end, ``flucut_report.html`` is then the dashboard of the
:doc:`live mode <live>`, which reloads itself in the browser every 30
seconds and shows the samples still to be processed as "not processed
yet" (the JSON summary has ``"running": true``, and ``summary.tsv`` and
``lengths.tsv`` only the samples already processed). The report and
the summary, their templates and plots come from Lumipore, with the
content of flucut.

The HTML report
---------------

The report is a single, self-contained file: the plots are inline SVG,
and it uses no external resources (and no JavaScript), so that it can
be opened offline in any browser, archived, and sent by e-mail. It
follows the light or dark theme of the browser.

Most figures of this chapter show the report of the simulated reads of
the benchmark (3,200 reads of the eight segments, of which 239 HA reads
internally primed, :doc:`performance`), trimmed with the references
(``-r``); the coverage figure is from a real Nanopore run, where it is
more telling. See ``docs/manual/figures/README.md`` for how they were
made.

With ``--report-segment`` (e.g. ``--report-segment HA``) the parts of
the report which list the references (the internal priming sites, the
complete amplicons and the coverage) show one segment only; the counts
of all reads (tiles, steps, primers, read ends by origin) are unchanged
With ``--segment``, which also limits the classification to the
references of one segment, the report shows that segment
(:doc:`cli_reference`).

**The run parameters.**
   At the right of the logo: the version of flucut, the command line,
   the primer scheme, the trimming policy, the references and the
   options of the trimming.

   .. figure:: figures/report_head.*
      :width: 100%

      The head of the report: the logo, the time of writing and the run
      parameters.

**Groups.**
   With several groups of input files (:doc:`output`), e.g. the barcode
   directories of a run, each group is a different sample and has its own
   section. The report starts with the throughput of the whole run (the
   tiles and the basic statistics) and a table of the groups, with their
   reads, the reads written, those with both primers and, with
   references, the internally primed read ends. Nothing else is pooled
   across groups: coverage, internal priming sites and amplicon lengths
   are shown for each group only. A column at the left of the page, which
   stays in view while the page scrolls, lists the groups: their names
   there, and in the table, link to their sections.

   .. figure:: figures/report_groups.*
      :width: 80%

      The groups of the test data ``testdata/fastq_pass`` (two barcode
      directories), without references.

Then, for each group (or for all reads, if there is only one):

**Sample: subtype and segments.**
   With references, what the located reads say about the sample itself,
   whatever is done with the trimmed reads afterwards. A first line
   gives the subtype (the H of the reference of its HA and the N of
   that of its NA) and the number of segments found, and notes when the
   segments are closest to references of different subtypes (a
   reassortant, or genes which differ little between subtypes). Then
   one line for each segment:

   - the reference its reads were located in, and the subtype named in
     its description;
   - the mean depth, the share of its positions covered by at least one
     read, and by at least 50 reads (what FLU-minion of IRMA needs to
     support a base);
   - the share of its depth that comes from internally primed reads,
     and from reads with an internal deletion: depth which does not
     represent the whole segment;
   - its complete amplicons, and the share of them shorter than the
     reference.

   A segment without located reads is shown as not found. Depth and
   coverage are those of the reads as flucut locates them (see Coverage
   by read origin below), not of a mapping.

**Reads at each step.**
   The reads read, those with the primers of both ends, of one end or
   of none, the primer dimers, the reads too short to be written, the
   reads written and those found on the minus strand, as tiles and as
   bars. Many reads without primers mean a scheme which does not match
   the primers of the run (``-s``), reads already trimmed, or a higher
   error rate than ``-e`` allows.

   Everywhere in the report, ``0.0%`` and ``100.0%`` mean exactly none
   and exactly all: a share which is not zero, but below 0.1%, is written
   ``<0.1%``, and one above 99.9%, but not all, ``>99.9%``.

   .. figure:: figures/report_tiles.*
      :width: 100%

      The tiles of all reads: 3,200 reads, 99.2% written, 88.4% with
      both primers, 231 internally primed read ends.

   .. figure:: figures/report_steps.*
      :width: 100%

      All reads: 99.2% of the reads were written; 88.4% had both
      primers, 11.6% only one (the simulation breaks some reads); the 24
      primer dimers were not written.

**Primers found.**
   The reads with the primer of each end, by primer of the scheme (e.g.
   the two forms of the Uni12 primer of the ``mbt`` scheme).

**Internal priming.**
   With references: the read ends of each end by origin (terminal,
   internal, not located), and the internal priming sites, by
   decreasing number of reads: the reference, the end, the position of
   the first base after the primer (1-based on the plus strand of the
   reference, as in the ``ml``/``mr`` tags of the reads) and the range
   of the positions of the site (the read ends of a site map a few bases
   apart, because of the sequencing errors: positions at most 5 bases
   apart are one site). Without references, a note: the sites are not
   known, and the length check below is the check of internal priming.

   .. figure:: figures/report_internal.*
      :width: 100%

      Internal priming in the simulated reads: 231 left ends were
      classified as internal; 224 of them are the site of the
      simulation, in the HA (KF918483.1), the first base after the
      primer at 1090.

**Complete amplicons.**
   For each reference (or for all reads, without references): the
   complete amplicons (both primers found), the most frequent length of
   the template between the primers, and the amplicons shorter than
   90% of the reference, with their most frequent length; then the
   histograms of the lengths, as in ``lengths.tsv``. Internal priming,
   and deletions, show as short amplicons.

   .. figure:: figures/report_complete.*
      :width: 100%

      Complete amplicons: 221 of the 367 HA amplicons are shorter than
      90% of the reference, most of them 670-679 bases long, from the
      internal site to the end of the segment.

   .. figure:: figures/report_lengths.*
      :width: 100%

      The histograms of the lengths of the complete amplicons, one for
      each reference and one for the amplicon whose read ends were not
      located; that of the HA (KF918483.1) has two peaks, the amplicons
      started at the internal site and the complete segment.

**Coverage by read origin.**
   With references: for each reference, a reconstructed coverage along
   its length, stacked by the origin of the reads covering each
   position. Each located read covers the span ``[start, start + length)``
   on the reference (from the mapped position of its probe, the
   ``ml``/``mr`` tags, and its kept length); the reads are coloured by
   where their amplicon started and on which strand:

   * **terminal +** and **terminal -**: started at a segment end (a
     genuine full-length amplicon), on the plus or the minus strand;
   * **internal +** and **internal -**: internally primed, on the plus
     or the minus strand;
   * **no 5' primer** and **5' probe not located** (grey, at the
     bottom): reads without a located 5' end. A read with no 5' primer
     is a fragment of unknown origin (a broken molecule, or a fragment
     of an Illumina library), placed by its 3' end and its length or,
     without a located 3' end, by its own sequence (see below).

   * **internal deletion** (purple, at the top): reads whose two ends
     map to the two ends of a segment, much farther apart than the read
     is long, e.g. defective interfering genomes. Such a read covers
     only its two flanks: the junction between them is located by
     probes from inside the read, so the band steps down where the
     deletions begin and up where they end. If the junction is not
     found, half of the read is drawn at each end.

   When only the 3' end of a read maps, its 5' origin is inferred from
   where its 5' end falls (the mapped 3' position less the read length):
   within the segment terminus it is counted as terminal, deeper as
   internal.

   A read of which neither end is located (no primer at all, or primers
   whose probes do not map) is still counted: a probe from the middle of
   its kept part is looked up in the k-mer index of the references, and
   the read covers the span around it given by its length. If the probe
   is not found, one at a quarter and one at three quarters of the read
   are tried. Only the probe is aligned, not the read, so the span is
   approximate (it ignores the indels of the rest of the read), but the
   depth is close to that of a read mapping: without this the reads
   without primers, most of an Illumina library, would be missing from
   the plot. The trimming does not depend on it. A read is left out only
   if no probe is found (a read which is not influenza, or one with an
   error in every k-mer of its probes).

   The same groups, with the strands merged, are the categories of the
   reads in the output: the ``ct`` tag, and the files written by
   ``--by-category`` (:ref:`reads-by-category`).

   A complete amplicon (both primers found) which could not be located
   in any of these ways is drawn over the whole segment of its length,
   as terminal, if the length tells which segment it is: within 80 bases
   of the length of a reference, and not as close to a reference of
   another segment (e.g. PB2 and PB1).

   Internal priming shows as a step up in the coloured (amber, blue)
   bands at the site, marked with the red line of the largest site of
   the reference. This is an overview built from flucut's classification,
   not a base-by-base mapping: it ignores indels. No other primer
   trimmer shows it.

   .. note::
      The coverage can only be reconstructed when the references are
      **complete segments, with the UTRs**. A read is located by the
      sequence just after its 5' primer, which in a full segment is the
      5' UTR; against coding-sequence (CDS) references, which lack the
      UTRs, only the internally primed reads map (their probe falls
      inside the ORF) and flucut warns (:doc:`troubleshooting`). Use
      ``scripts/fetch_references.py`` for complete segments.

   .. figure:: figures/report_coverage.*
      :width: 100%

      The coverage by read origin of a real Nanopore run of an H5N1
      isolate (SRA ``SRR36760543``, against complete H5N1 segment
      references): every segment is covered by full-length amplicons
      (green), and the HA (segment 4) shows the internal-priming step
      at the site (the red line) with the internally primed reads (amber,
      blue) rising far above the terminal baseline — internal artefacts
      (:doc:`internal_priming`). PB1 shows a smaller internal step.

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
  group, with ``reads`` (the reads at each step), ``primers``,
  ``ends`` and ``internal_sites`` (with references), ``short_complete``
  and ``complete_amplicons`` (with the histograms of the lengths).

For example, the internal priming sites of each group with ``jq``::

   jq -r '.groups[] | .name as $g | .internal_sites[]? |
          [$g, .reference, .end, .position, .reads] | @tsv' \
       out/flucut_summary.json

which gives for the simulated reads::

   reads  KF918483.1   left  1090  224
   reads  KF918483.1   left  1158  6
   reads  NC_002023.1  left  1595  1
