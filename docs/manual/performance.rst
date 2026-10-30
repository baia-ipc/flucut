Performance and accuracy
========================

The numbers of this chapter come from the benchmark of the repository
(``benchmark/``): what was run, on which data, with which versions and
commands, is in ``benchmark/RESULTS.md``; the scripts are in
``benchmark/``. The reads are simulated, so that the origin of every
base is known: amplicons of a genome made of the segments of
A/Puerto Rico/8/1934 (H1N1) and the HA of A/chicken/Cambodia/X0815302/
2013 (H5N1), which carries the internal priming site of the I348K
artefact, with Nanopore-like errors (``benchmark/README.md``). No real
run has been measured yet.

The consensus
-------------

On 3200 simulated reads (400 per segment, 60% of the HA reads started
at the internal site), assembled by IRMA 1.4.0 (FLU-minion):

- the untrimmed reads give the artefact: AAA (lysine) at the I348K
  codon of the majority consensus; the reads trimmed by flucut give
  ATA, the truth;
- with the reference elongation of IRMA, the consensus of the reads
  trimmed by flucut has no error, and ends 12 and 13 bases short of the
  segments (the primers); with ``--trim terminus`` it has complete ends,
  with the primer's bases there (:doc:`irma`);
- at 10, 20 and 40 reads per segment (three simulations each), the
  untrimmed reads give the artefact in 5 of 9 runs, the reads trimmed by
  flucut in none.

The reads
---------

Bases of the 3200 reads by their true origin, after the trimming (flucut
of the current version; the consensus above was made with the reads of
the previous one, which removed 882 template bases, kept 1852 primer
bases and 7557 flank bases):

=================== ========= ========= ======== ===========
tool                template  primer    flank    internal
                    removed   kept      kept     site kept
=================== ========= ========= ======== ===========
untrimmed           0         102554    127032   2698
flucut              680       1819      7542     23
cutadapt            102408    3037      8784     25
=================== ========= ========= ======== ===========

cutadapt (``-g`` and ``-a``, not anchored, as the reads have flanking
bases) removes everything before a primer-like sequence found anywhere
in a read: with 2 errors, Uni12 also matches inside the segments, e.g. at
the I348K site, and complete HA reads are cut there. Its consensus is
nonetheless the same as that of flucut in these simulations, also at low
coverage. The comparison with the other tools (IRMA-core trimmer,
BBDuk, Porechop_ABI, dorado, Barbell, samtools ampliconclip) is in
``benchmark/RESULTS.md``.

Reads with barcodes
-------------------

On 3200 reads simulated with the native barcoding adapters of Oxford
Nanopore (NB01-NB04) outside the primers:

- the untrimmed reads gave the right codon at the I348K site, but only
  by 158 reads against 152: the aligner of IRMA clips more of the start
  of the internally primed reads after an adapter. The artefact thus
  depends on the library and on chance; after flucut, 157 reads against
  1;
- with the elongation, the consensus of the reads trimmed by flucut has
  no error, as with cutadapt, dorado, samtools ampliconclip and Barbell;
- without the adapter search (``--adapters none``), flucut removed 254
  template bases and left 2221 primer bases and 28539 bases of adapters,
  at the read ends where it found no primer; with it (the default), 267,
  1019 and 698 (Barbell: 1883, 445 and 4375; dorado: 3714, 6452 and 360).
  The consensus is the same. The reads without barcodes are trimmed as
  without the search.

Speed and memory
----------------

On 32000 reads (53 million bases; ``benchmark/RESULTS.md``), with 4
threads on a 4-core machine: flucut 1.6 s (61 MB; 3.8 s and 27 MB with
one thread; with 8 references, 1.8 s and 78 MB), cutadapt 3.7 s (41
MB), BBDuk 4.0 s (954 MB), dorado trim 8.1 s (479 MB), minimap2 with
samtools ampliconclip 15 s (125 MB), Porechop_ABI 34 s (253 MB); the
IRMA-core trimmer, single-threaded, 7.5 s (11 MB); Barbell, on barcoded
reads, 2.4 s (44 MB; it writes uncompressed FASTQ).

A large part of the time of flucut is the compression of its output
(gzip, one member for each batch of 1000 reads). It is compressed with
libdeflate: with zlib, at the same level (6), flucut took 2.4 s with 4
threads and 6.6 s with one, for an output 4% larger. The level is set
by ``--compress-level`` (on the same reads; decompressing the input
alone takes 0.8 s):

.. list-table::
   :header-rows: 1

   * - compression
     - 1 thread, wall s
     - 1 thread, CPU s
     - 4 threads, wall s
     - 4 threads, CPU s
     - output MB
   * - zlib, level 6 (before)
     - 6.58
     - 6.41
     - 2.43
     - 7.98
     - 41.0
   * - libdeflate, level 1
     - 3.28
     - 3.14
     - 1.28
     - 3.76
     - 49.1
   * - libdeflate, level 4
     - 3.64
     - 3.40
     - 1.50
     - 4.80
     - 42.2
   * - libdeflate, level 6 (default)
     - 3.90
     - 3.70
     - 1.63
     - 4.66
     - 39.3
   * - libdeflate, level 9
     - 5.56
     - 5.34
     - 2.26
     - 7.10
     - 37.2

The earlier measures
--------------------

flucut trimmed the 3200 reads (5.3 million bases) in 0.7 seconds with
one thread. The classification against references, through the k-mer
index, adds little: 0.8 s instead of 0.7 s with 8 references, 1.1 s with
80 (measured before the search of the truncated primers; with the
exhaustive search, 2.0 and 14-16 s).

With threads: the same reads ten times (32,000 reads), with the 8
references, on a computer with 4 cores (three runs each): 6.7-6.9 s
with ``-j 1``, 2.1-2.5 s with ``-j 4`` (3.0-3.2 s with the threads of
the previous version, started for each batch); the output files are
identical.
