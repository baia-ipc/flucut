Quick start
===========

The repository contains a small test set in ``testdata``: twelve reads
simulated from a real HA record (NCBI ``NC_004908.1``), each one a
typical situation (see ``testdata/README.md``). From the root of the
repository, after ``make``:

::

   bin/flucut -o out testdata/reads.fastq

In a terminal, flucut shows its banner, then:

::

   SUCCESS | out/reads.fastq.gz: 11 of 12 reads written (10 with both
     primers, 1 with one, 1 with none; 1 too short, of which 1 dimers)
   SUCCESS | out/reads.fastq.gz: complete amplicons, lengths between the
     primers (peaks): 1700-1709 (5) 1000-1009 (2) 20-29 (1) 900-909 (1)

(the lines are wrapped here). The output directory contains:

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

::

   SUCCESS | out2/reads.fastq.gz: ends started by internal priming: 2 left,
     2 right; not located in the references: 1 left, 1 right
   SUCCESS | out2/reads.fastq.gz: NC_004908.1: 4 of 9 complete amplicons
     shorter than 90% of the reference (1714 bases), most at 1000-1009:
     internal priming or deletions

and the read simulated as started inside the segment, at position 712,
is tagged as such:

::

   @internal12_plus simulated  fs:Z:+  fl:Z:uni12r  ol:Z:internal
     ml:Z:NC_004908.1:712:+  fr:Z:mbtuni-13  or:Z:terminal
     mr:Z:NC_004908.1:1714:-  tl:i:12  tr:i:23

To trim a real run, give the reads of each sample, as files or as
directories (one output per input):

::

   flucut -o trimmed barcode01/ barcode02/ barcode03/

and give the outputs to IRMA (:doc:`irma`).
