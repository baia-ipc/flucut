Quick start
===========

The repository contains a small test set in ``testdata``: twelve reads
simulated from a real HA record (NCBI ``NC_004908.1``), each one a
typical situation (see ``testdata/README.md``). From the root of the
repository, after ``make``:

::

   bin/flucut -o out testdata/reads.fastq

flucut writes:

::

   SUCCESS | out/reads.fastq.gz: 11 of 12 reads written (10 with both
     primers, 1 with one, 1 with none; 1 too short, of which 1 dimers)
   SUCCESS | Reports written: out/flucut_report.html,
     out/flucut_summary.json

(the lines are wrapped here, without the time stamps). The output
directory contains:

``reads.fastq.gz``
   the trimmed reads, ready for IRMA;
``summary.tsv``
   the counts of the reads;
``lengths.tsv``
   the lengths of the complete amplicons.

The first read, an amplicon with both primers, became:

::

   @full_plus simulated  fs:Z:+  fl:Z:uni12r  fr:Z:mbtuni-13  tl:i:12  tr:i:23
   GGAATTACTTAACTAGCAAAATGGAAACAATATCACTAATAACTATACTACTAG...

The tags tell that the read is of the plus strand of the amplicon
(``fs``), that the primer ``uni12r`` was found at its left end and
``mbtuni-13`` at its right end, and that 12 and 23 bases were removed
at its start and end (:doc:`output`). The primer dimer (two primers and
nothing between them) was not written.

With the reference of the segment, flucut also tells where each read
started:

::

   bin/flucut -o out2 -r testdata/ref-ha.fasta testdata/reads.fastq

It reports the ends started by internal priming and the complete
amplicons much shorter than the reference (the figure shows this run in
a terminal, where flucut shows its banner first), and the read simulated
as started inside the segment, at position 712, is tagged as such:

::

   @internal12_plus simulated  fs:Z:+  fl:Z:uni12r  ol:Z:internal
     ml:Z:NC_004908.1:712:+  fr:Z:mbtuni-13  or:Z:terminal
     mr:Z:NC_004908.1:1714:-  tl:i:12  tr:i:23

.. figure:: figures/terminal_flucut.*
   :width: 100%

   ``flucut -o out -r testdata/ref-ha.fasta testdata/reads.fastq`` in a
   terminal of 100 columns (``terminal.js``).

To trim a real run, give the reads of each sample, as files or as
directories (one output per input):

::

   flucut -o trimmed barcode01/ barcode02/ barcode03/

or the ``fastq_pass`` directory of the run, which gives one output for
each of its barcode directories:

::

   flucut -o trimmed run/fastq_pass/

and give the outputs to IRMA (:doc:`irma`).
