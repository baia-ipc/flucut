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

HA2 begins directly after the cleavage site with the fusion peptide,
which every haemagglutinin has and which lies outside the
low-complexity stretch. It reads ``GLFGAIAG`` in most of them,
``GIFGAIAG`` in the H3 of horses and pigs, ``GLLGAIAG`` in some others
and ``GFFGAIAG`` in influenza B, so flucut looks for the pattern
``G..GAIAG``, leaving the two residues which vary free. It occurs
exactly once in each of the 84 HA references of the reference sets
flucut ships, and never more than once in 4484 HA records of every
subtype H1 to H16 and of influenza B downloaded from GenBank, so the
anchor cannot be found in the wrong place in a read. In 39 of those
records it is not found at all, always because the peptide is absent
and never because its placement is ambiguous.

The pattern is deliberately no looser, and its sixth position in
particular is fixed. Twenty-one of those records carry a lysine there
instead of the isoleucine, spelled ``AAA`` in every one, across eight
subtypes: that is what internal priming by the Uni12 primer writes over
the start of HA2 (I351K, residue 6 of HA2, see
:doc:`internal_priming`), deposited as if it were sequence. A pattern
which accepted it would read a cleavage site out of records whose
fusion peptide is a primer.

flucut takes those 24 bases of the HA reference as an *anchor*, finds the anchor inside each read (in
both orientations, with the error allowance of the probes), and reads
the 24 bases immediately before the match, as the read carries them.

Because the anchor sits directly beside the site, an insertion or a
deletion anywhere else in the read cannot shift what is read out of
frame. Nothing is mapped and nothing is assembled: one search of a
24-base pattern in each read is enough, so the site is reported while
the run is still sequencing.

flucut also reports the fusion peptide itself, as the reads carry it
where the anchor matched. It costs nothing, since those bases are
already in hand, and it answers two questions the site alone cannot: the
peptide is conserved, so reads which agree on it but not on the site
have a genuinely hard site, while reads agreeing on neither have an
anchor matching reads it should not; and a sample whose own fusion
peptide differs from its reference's belongs to a lineage that reference
does not represent.

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

Two things are measured from the same reads, and both must hold.

First, *do the reads agree at the site?* The share of them carrying the
majority residues is the **support**.

Second, *do they agree where nothing is difficult?* Beside the anchor,
on the HA2 side, lie two windows of the same width as the site, taken
from the same reads, where HA2 is ordinary coding sequence. The share
agreeing there is the **control support**. It tells a site the reads
cannot resolve from reads that agree nowhere: a sample high in the
control and low at the site has a genuinely hard site, while a sample
low in both has reads too noisy, or holding more than one virus, for a
claim about any position.

.. list-table::
   :header-rows: 1

   * - measured
     - needed to be firm
   * - reads over the site
     - 20
   * - support at the site
     - 0.45
   * - support in the control windows
     - 0.60

The verdict is one of

``firm``
   enough reads, enough of them agree at the site, they agree beside it,
   and they read like the reference of the sample.

``confirm``
   one of those fails. The residues are still reported, with the reason.
   Confirm before reporting a pathotype.

``too few reads``
   fewer than 20 reads reach the site; nothing is said about it.

Neither threshold names a platform or a chemistry. The chemistry is
still read from the basecalling model in the read titles
(``dna_r10.4.1_…``, ``dna_r9.4.1_…``) and reported, because it is worth
knowing, but it decides nothing: on 104 samples the control windows
separated the chemistries by themselves, every R9 sample falling below
the floor and the R10 samples sitting well above it. A bad flow cell, or
a chemistry that does not exist yet, is then judged on the reads in
front of it rather than on a label.

The thresholds are the numbers of one laboratory's runs and are a
starting point: a laboratory whose runs differ should look at what its
own clean samples reach before trusting the line between ``firm`` and
``confirm``. Everything reported beside the verdict — the residues, the
reads behind them, both supports and the other sites found — does not
depend on them.

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

How many sites the reads hold is the number of **supported sites**: the
residues carried by at least 10% of the reads (and at least 2 reads).
The number of distinct residues seen in the reads (``distinct_sites`` in
the JSON summary) and of distinct spellings (``distinct_motifs``) is not
that: a read with an error has a reading of its own, so both count the
errors of single reads and are large in nearly every sample. Measured on
the 104 Nanopore samples below, 10% leaves one site in 89 samples and
two or three in 13, most of them holding a second haemagglutinin; at 5%
or less the errors that recur in the same sample, one residue wrong in
5 to 13% of the reads, are counted too.

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

How this compares with other tools
----------------------------------

Several tools report the cleavage site, and all of them read it from an
assembled consensus: nf-flu counts basic residues in the consensus
annotated by VADR and labels the sample LPAI or HPAI; the Nextstrain
avian-flu workflow finds ``GLFG`` in the HA alignment and translates the
bases before it; the Nextclade A(H5) datasets flag a fixed window of HA
residues as polybasic when it holds four or more R or K. None of them
reports how many reads support the call. IRMA, GenoFLU, LABEL, INSaFLU
and FluMut do not report the site at all.

The idea of finding the site by the fusion peptide that follows it is
not new — the Nextstrain script does the same — and what flucut adds is
to do it one read at a time, before an assembly exists, and to say how
many reads agree.

The nearest work takes the opposite route. DeGenRESOLVE repairs the
consensus itself, adjudicating the indels a variant caller proposes by
depth, by the fraction of reads behind them, and by a frame test. Its
stated limit is where this check begins: its result is a dominant
consensus, not a record of what the reads disagreed about. The two are
complementary.

See ``docs/cleavage-site.md`` for the comparison in full, with
references.

What it is not
--------------

It does not assign a pathotype. A multibasic cleavage site is the
sequence determinant of high pathogenicity, not a diagnosis: the
pathotype is decided by the designated authorities, on more than the HA
segment and usually with a test in birds.

It does not correct a consensus. A majority over reads which gain and
lose bases in this stretch is itself biased towards losing them, so a
site supported by few reads is a reason to confirm it, not an answer.

It reads eight codons. That covers the basic run of every reference
flucut ships, the longest of which is six residues, but an insert
longer than that, as some highly pathogenic H7 carry, is shown only in
part; the site still reads as multibasic, so only the residues printed
are cut short.

It does not assign a clade. The site and the clade are read in different
ways and at different moments: the site from single reads beside the
fusion peptide, before assembly; the clade from a near-complete
haemagglutinin, after it. Assemble and then run a clade classifier on
the consensus — Nextclade with an avian A(H5) dataset, or LABEL, which
IRMA ships. :doc:`markers` says more about why the boundary is drawn
there.

Example
-------

.. code-block:: text

   HA cleavage site

   QRERRRKR           ✓ firm   · multibasic

   reads it was read from        6,977 of 9,588 searched
   reads with these residues     4,524 (64.8%)
   sites in at least 10% of the reads  1
   reads with a multibasic site  92.3%
   reads agreeing beside the site 87.7% of 13,410
   chemistry                     R10
   fusion peptide in the reads   GLFGAIAG (99.1% of 6,977)
   reference: site and peptide   LREKRRKR, GLFGAIAG

A sample holding the reads of a second haemagglutinin instead reads

.. code-block:: text

   HA cleavage site

   NVPSRSSR           ! confirm · monobasic

   reads it was read from        65 of 950 searched
   reads with these residues     41 (63.1%)
   sites in at least 10% of the reads  2
   reads with a multibasic site  30.8%
   reads agreeing beside the site 61.7% of 128
   chemistry                     R10
   fusion peptide in the reads   GLFGAIAG (96.9% of 65)
   reference: site and peptide   LREKRRKR, GLFGAIAG

   the 5 most frequent sites in these reads
     NVPSRSSR  41  63.1%  monobasic
     QKERRRKR   9  13.8%  multibasic
     LRERRRKR   3   4.6%  multibasic
     …

   Most of the reads do not read like the reference of this sample
   (LREKRRKR): it holds the reads of more than one haemagglutinin, and
   the site above is not necessarily the one of its own virus.
