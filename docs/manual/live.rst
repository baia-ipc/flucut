Live mode
=========

With ``--live``, flucut follows a MinKNOW run while it is sequencing: the
files of reads are trimmed as they are written, and the outputs, the
tables and the reports are updated after each file, so that the samples
can be assembled, or at least looked at, before the end of the run.

::

   flucut --live -o out /data/run1/fastq_pass

The inputs must be directories: typically the ``fastq_pass`` directory of
the run, whose ``barcode`` subdirectories are the samples (also those
created during the run), or a directory of files of one sample.

How the files are taken
-----------------------

The input directories are scanned every ``--scan-interval`` seconds
(default 10). A new file is trimmed when it is ready: not modified for
``--min-age`` (default 60 s; e.g. ``30s``, ``2m``) and complete (a
gzipped file whose stream is not truncated, a FASTQ file whose last line
is complete, a BAM file ending with its end-of-file marker). Its reads
are appended to the output of its sample (``<sample>.fastq.gz``: gzip
files can be concatenated), and one line per file reports them, e.g.:

::

   SUCCESS | barcode01: FBA12345_pass_barcode01_0.fastq.gz: 3912 of 4000
     reads written (3801 with both primers, ...)

A file which changes after it was trimmed is reported (``modified after
being processed``): its new reads are not trimmed; ``--min-age`` should
then be increased.

``-e auto`` (the default) is decided at the first file ready, from its
first 1000 reads, as in a normal run (:doc:`cli_reference`). So is the
scheme with ``-s auto`` (the default, :doc:`primers`), and the platform
with ``-p auto`` (the default, :ref:`platform-check`), before anything
else: the settings, the counts and the reports depend on them, so the
mode waits for that file before writing the first report. A platform
given is checked against the first file processed.

The reports
-----------

After each file, ``summary.tsv``, ``lengths.tsv``, the JSON summary and
the dashboard, ``flucut_dashboard.html``, are written again (each to a
temporary file, then renamed, so that a reader never sees an incomplete
one). flucut says where the dashboard is when it starts watching.

The **dashboard** (:doc:`reports`) reloads itself in the browser (the
JSON summary has ``"running": true``). The reload interval adapts to the
run: it is halved when a scan has processed a new file and doubled when
a scan finds none, between 1 and 30 seconds, so the page refreshes
quickly while reads are coming in and rarely once the run is quiet. As
each barcode directory is a different sample, the dashboard does not
pool them: it shows the overview of the run, the HA in all samples, or
one sample. The choice is kept in the page address, so it survives the
automatic reload. While reads are still being added, the count of
complete segments of a sample is grey until all are complete: a
segment not yet at the depth of 50 is not a failure during the run.

.. figure:: figures/dashboard.*
   :width: 100%

   The dashboard of a run of four barcodes (simulated reads).

At the end, the dashboard stays in its last state and the **report**
(:doc:`reports`) is written, with a section for each barcode; the counts
of each sample are reported as in a normal run, with the warnings of the
lengths.

The end
-------

The live mode ends:

- at the end of the run: when MinKNOW writes its final summary
  (``final_summary_*.txt``) in an input directory or in its parent (the
  run directory of ``fastq_pass``); the files present are then trimmed
  without waiting for their minimum age;
- after ``--timeout`` without new files (e.g. ``30m``, ``2h``; a number
  is in minutes; default: none);
- with Ctrl-C (or SIGTERM): the file being trimmed is completed, then
  the final reports are written; a second Ctrl-C stops flucut at once.

Resuming
--------

The list of the files trimmed for each sample and its counts are saved
in the output directory (``.flucut-state/<sample>.flucut-files``,
``.flucut-state/<sample>.flucut-counts``) after each file. With
``--resume``, a live run started again with the same output directory
skips the files already trimmed, and continues the outputs and the counts: e.g. after a
restart of the computer. Without ``--resume``, the outputs of each
sample are written again from the beginning.

The outputs of a live run are the same as those of a normal run on the
files of the end of the run (``tests/test-live.sh`` checks it, also
after an interruption and a resumption).
