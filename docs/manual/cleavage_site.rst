The HA cleavage site
====================

What it is
----------

The haemagglutinin is made as one protein and cut in two, HA1 and HA2,
before the virus can fuse with a cell. The residues at the cut are the
*cleavage site*. In the highly pathogenic avian H5 and H7 viruses it is
a run of basic residues (arginine R, lysine K), for example
``PQRERRRKR``; in the low pathogenic ones a single one, for example
``PQRETR``. The difference is what distinguishes the two pathotypes, so
a consensus which is wrong there describes the wrong virus.

The run of basic residues is encoded by an A/G-rich stretch of the
segment. Low-complexity sequence of that kind is where sequencing reads
most often gain or lose a base, and an assembler working from few reads
writes such an error into the consensus. The reading frames of HA1 and
HA2 then differ and the site reads as nonsense, or, with three bases too
many, the run of basic residues is one residue longer than it should be
while the frame still looks intact.

flucut reports what the **reads** carry at the site, and how many of
them agree. It is not an assembly and does not replace one: it says what
the data supports before any consensus is built, and whether the
consensus that follows can be trusted there.

How flucut reads it
-------------------

HA2 begins directly after the cleavage site with the fusion peptide
``GLFGAIAG``, which is the same in every influenza A haemagglutinin and
lies outside the low-complexity stretch. flucut takes those 24 bases of
the HA reference as an *anchor*, finds the anchor inside each read (in
both orientations, with the error allowance of the probes), and reads
the 24 bases immediately before the match, as the read carries them.

Because the anchor sits directly beside the site, an insertion or a
deletion anywhere else in the read cannot shift what is read out of
frame. Nothing is mapped and nothing is assembled: one search of a
24-base pattern in each read is enough, so the site is reported while
the run is still sequencing.

The reads are counted by the **residues** they carry, not by the bases.
The residues decide the pathotype, and the reads agree on them far more
often than on the spelling, which differs by a base wherever a read has
an error.

Which reads it reads
--------------------

What reaches the site differs between the platforms:

* with ``-p nanopore`` (the default), a read covers the whole amplicon,
  and it reaches the site only when its amplicon started at a segment
  end. The internally primed reads begin just after the site, because
  the Uni12-like site lies a few bases after it, and the genomes with an
  internal deletion have removed that part of the segment. flucut
  therefore reads the site from the **full templates** located on the HA
  reference of the sample.
* with ``-p illumina``, a read is a fragment of an amplicon. A fragment
  which covers the site carries no primer, so it has no category and is
  not located at all. flucut therefore searches **every read**, and the
  anchor itself says which ones cover the site.

How far to trust it
-------------------

Two numbers decide: how many reads gave the site, and what share of them
agree on the same residues, the *support*. How much support is enough
depends on how often the reads gain or lose a base in this stretch,
which differs between the kinds of data:

.. list-table::
   :header-rows: 1

   * - data
     - support needed to be firm
   * - Illumina
     - 0.60
   * - Nanopore R10
     - 0.45
   * - Nanopore R9
     - never firm on its own

The chemistry is read from the basecalling model in the titles of the
reads (``dna_r10.4.1_…``, ``dna_r9.4.1_…``). When the titles do not say
it, flucut uses the R10 number and reports the chemistry as unknown.

The verdict is one of

``firm``
   enough reads, enough of them agree, and they read like the reference
   of the sample.

``confirm``
   the residues reported are the majority, but too few of the reads
   agree for them to stand on their own, or the reads do not read like
   the reference of the sample. Confirm before reporting a pathotype.

``too few reads``
   fewer than 20 reads reach the site; nothing is said about it.

Two numbers beside the majority are often more useful than it:

* the share of the reads whose **own** site is multibasic, whatever its
  spelling. The reads can disagree on the exact residues and still agree
  that the site is multibasic, which is what the pathotype turns on;
* the **other sites** found in the same reads. The fusion peptide is the
  same in every subtype, so the anchor finds the site of every
  haemagglutinin in the sample. A barcode holding the reads of more than
  one of them shows two sites with substantial shares; when the majority
  does not read like the HA reference chosen for that sample, flucut
  says so and never calls the result firm.

Where it is reported
--------------------

In the HTML report and the dashboard, as a card in each sample and as a
column of the table of the HA in all samples; in ``summary.tsv`` as the
columns ``cleavage_site``, ``cleavage_reads``, ``cleavage_support``,
``cleavage_multibasic`` and ``cleavage_verdict``; and in
``flucut_summary.json`` as the ``cleavage_site`` object of each section.

It needs references, so it is not reported with ``-r off``, and it needs
an HA reference which contains the fusion peptide. A reference which is
a coding sequence is enough, as the site lies inside the reading frame.

What it is not
--------------

It does not assign a pathotype. A multibasic cleavage site is the
sequence determinant of high pathogenicity, not a diagnosis: the
pathotype is decided by the designated authorities, on more than the HA
segment and usually with a test in birds.

It does not correct a consensus. A majority over reads which gain and
lose bases in this stretch is itself biased towards losing them, so a
site supported by few reads is a reason to confirm it, not an answer.

Example
-------

.. code-block:: text

   HA cleavage site

   QRERRRKR           ✓ firm   · multibasic

   reads it was read from        6,977 of 9,588 searched
   reads with these residues     4,524 (64.8%)
   reads with a multibasic site  92.3%
   chemistry                     R10
   site of the reference         LREKRRKR

A sample holding the reads of a second haemagglutinin instead reads

.. code-block:: text

   HA cleavage site

   NVPSRSSR           ! confirm · monobasic

   reads it was read from        65 of 950 searched
   reads with these residues     41 (63.1%)
   reads with a multibasic site  30.8%
   chemistry                     R10
   site of the reference         LREKRRKR

   the 5 sites in these reads
     NVPSRSSR  41  63.1%  monobasic
     QKERRRKR   9  13.8%  multibasic
     LRERRRKR   3   4.6%  multibasic
     …

   Most of the reads do not read like the reference of this sample
   (LREKRRKR): it holds the reads of more than one haemagglutinin, and
   the site above is not necessarily the one of its own virus.
