Markers
=======

What they are
-------------

A few positions of the influenza proteins are worth knowing before
anything is assembled: the residues that decide whether a virus binds
the receptors of birds or of humans, whether its polymerase works in a
mammalian host, and whether it still responds to the antivirals. flucut
reads those positions straight from the reads, as it reads the
haemagglutinin cleavage site, and reports what the reads carry at each
of them.

The report is a distribution, not a call. For each marker you get the
residue, the number of reads it was read from, and the share of those
reads carrying each value. A residue present in only a part of the
sample is therefore visible as such: a consensus would round it away,
and with it the first sign that a part of the virus population has
changed.

What flucut reads
-----------------

.. list-table::
   :header-rows: 1
   :widths: 14 8 16 24 38

   * - marker
     - protein
     - positions
     - numbering
     - what the residue means
   * - ``PB2-627``
     - PB2
     - 627
     - PB2, from the initiator methionine
     - lysine replicates in mammalian cells at their temperature;
       glutamate is the avian residue
   * - ``PB2-701``
     - PB2
     - 701
     - PB2, from the initiator methionine
     - asparagine adapts the polymerase to mammals, as 627K does and
       without it
   * - ``PA-38``
     - PA
     - 38
     - PA, from the initiator methionine
     - threonine lowers the susceptibility to baloxavir
   * - ``NA-275``
     - NA
     - 275
     - N1 numbering (274 in N2 numbering)
     - tyrosine lowers the susceptibility to oseltamivir in the N1
       neuraminidases
   * - ``HA-226-228``
     - HA
     - 226, 227, 228
     - H3 numbering of the mature HA1
     - leucine at 226 with serine at 228 binds the receptors of the
       human airway; glutamine with glycine binds the avian ones
   * - ``M2-31``
     - M2
     - 31
     - M2, from the initiator methionine
     - asparagine resists amantadine and rimantadine, and is the
       residue of most viruses now circulating
   * - ``NA-stalk``
     - NA
     - the stalk, residues 36 to 90
     - bases from the start of NA to its conserved HxNxT peptide
     - a shorter stretch is a deletion in the stalk, which
       poultry-adapted viruses carry
   * - ``NS1-tail``
     - NS1
     - the C-terminus
     - the residues before the stop codon
     - the last four residues are the PDZ-binding motif, and a stop
       before them is the truncated C-terminus of the viruses
       circulating in humans
   * - ``PB1-F2``
     - PB1-F2
     - the whole reading frame
     - codons of PB1-F2
     - a frame of 87 to 90 codons encodes the whole protein; a stop
       codon early in it leaves a fragment

These are markers of influenza A. No anchor of the table is found in the
influenza B references, so an influenza B sample reports no marker
rather than a wrong one.

What is not in the table, and will not be: clades, genotypes and
lineages. They are defined by nomenclatures which are revised, and a
table compiled into a program cannot follow them; a program which named
a clade would be wrong a year later without anything in it having
changed. A marker is listed only when its evidence is settled and when
a short conserved peptide lies directly beside it, so that the position
can be found in the sequence itself.

How flucut finds them
---------------------

Not by counting residues from the start of a protein. The position of a
residue depends on what it is counted on: the haemagglutinin is numbered
on its precursor, on the mature protein or on the structure of H3, the
neuraminidase of N1 and of N2 differ by one at the oseltamivir position,
and a virus whose stalk is shorter has every residue after it shifted.
In the reference set flucut ships, the residue equivalent to H275 of N1
sits at position 272, 274 or 275 of the reference depending on the
length of its stalk.

Instead, each position is found by its **anchor**: the bases of the
reference which encode a short conserved peptide directly beside it. The
anchor is searched in each read, in both orientations and with errors
allowed, and the residues beside the match are read as the read carries
them. Because the anchor is directly beside the position, an insertion
or a deletion elsewhere in the read cannot shift what is read out of
frame, and a virus whose numbering is shifted still gives the right
residue.

An anchor is used for a reference only when the peptide occurs exactly
once in it and the bases of the anchor are found nowhere else in it
within the errors a read is searched with. Both are checked when a run
starts; a reference which fails either is named at ``--debug`` and
simply gives no marker.

Finding the anchor is not enough: the bases beside it have to be the
bases of that position. A marker is therefore read only from a stretch
of a read which follows its reference — within the part flucut keeps,
which is the part it judges to be the template, clear of the junction
of an internal deletion and of any primer, and at a place the read
actually covers. The bases flucut removes are the whole primer at an
internally primed end, and reading over them reads the primer: the
conserved terminus ``AGCAAAAGCAGG`` ends in ``AAA``, a lysine, so a
marker whose window is the codon before its anchor would read K on
every read primed just before it.

How much to trust a value
-------------------------

Each marker comes with three numbers and a verdict.

**reads**
   How many reads the value was read from. The markers in the middle of
   a segment are read from fewer reads than those near its ends,
   especially in a sample full of defective genomes, which keep the two
   ends and lose the middle.

**support**
   The share of those reads carrying the value reported.

**anchor support**
   The share of the same reads agreeing on the residues where the anchor
   matched. The anchor is conserved, so this says how much the reads of
   this sample agree where nothing is difficult — the floor below which
   a value carried by few reads cannot be told from the errors of the
   reads.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - verdict
     - meaning
   * - ``firm``
     - at least 60% of the reads agree on the value, and they agree at
       the anchor
   * - ``confirm``
     - the value is the majority, but too few reads agree for it to
       stand alone
   * - ``unresolved``
     - fewer than 35% of the reads agree on any value: none is reported
   * - ``too few reads``
     - fewer than 20 reads

A value carried by at least a tenth of the reads and by at least 3
reads counts as one the sample holds. Where there is more than one, the
sample holds a mixture, and the report shows each value with its share.
Below a tenth a value is not told from the disagreement the reads
already show at the anchor, where nothing is difficult: on 102 Nanopore
samples that disagreement is a median of 14% of the reads on R10 and
39% on R9.

This is where reading the reads differs from reading a consensus. On
those samples 45 residues are carried by at least a tenth of the reads
and are not the one in the consensus of the same reads — in four
samples the amantadine residue M2 31N, in 41% of 22,645 reads in one of
them, where the consensus shows only serine.

``PB1-F2`` and ``NA-stalk`` are read over far more than one codon, so
one indel anywhere in the window changes what is read. On Nanopore
reads they are often ``unresolved``, or read from too few reads; on
Illumina reads, whose indel rate is far lower, they are not. They are
reported that way rather than left out.

The verdict is about the value, not about the virus: ``firm`` means the
reads agree, not that the sample holds one virus.

Reading the report
------------------

In the HTML report and in the dashboard, each sample has a **Markers**
card: one row per marker, with the value the reads carry, the value of
the reference, the share of the reads agreeing and the verdict. A badge
marks a value worth a look. Where the reads carry more than one value at
a marker, the values and their shares are listed under *the markers at
which the reads carry more than one value*, and the table *what each
marker is* gives the positions, the numbering convention and the
meaning.

The table of a segment across the samples has a **markers** column, so
that a run can be scanned for the samples whose residues differ.

In ``summary.tsv`` there is a column per marker, named after it
(``pb2_627``, ``ha_226_228``, ...), holding the value only where the
reads support one, and a column ``markers_attention`` naming the markers
of that sample which deserve a look. The full detail, with the share of
the reads carrying each value, is in the ``markers`` array of
``flucut_summary.json``.

A sample whose reads carry a residue worth a look is flagged with it, so
that it is seen in the run overview and not only in the sample. The
residues flagged are 627K, 701N, 38T, 275Y, 226L, 228S and a stalk
shorter than the one of the reference. M2 31N is reported but not
flagged: nearly every virus now circulating carries it, so the remark
would be on every sample.

What this is not
----------------

A marker is a residue read from the reads, nothing more.

* It is not a phenotype. A residue associated with resistance in the
  literature is reported as that residue, not as resistance; a
  susceptibility result needs a susceptibility assay.
* It is not a clade, a genotype or a lineage.
* It does not replace a consensus: a marker is a handful of positions
  read early, and the rest of the segment still has to be assembled.

Why there is no clade, and what to run instead
----------------------------------------------

If you see the marker residues and the cleavage site in the report, the
next question is usually where the clade is. It is not there, and it is
not coming.

A marker is a residue; a clade is a placement on a tree. The residue can
be read from one read, because a short conserved peptide sits directly
beside it and says where it is. A clade has no such peptide: it is
decided by the pattern of a whole haemagglutinin against the
haemagglutinins of every other virus, so it needs a near-complete
segment, which exists only once the reads have been assembled. The
anchor flucut uses needs a fixed position in conserved context, which is
the opposite of a site that is informative because it varies.

There is a practical reason as well. Clade nomenclatures are revised,
and a table built into a program cannot follow them: it would name a
clade correctly this year and wrongly the next without anything in the
program having changed. The tools that do assign clades are maintained
with their nomenclature, which is what makes them right.

So assemble first and then run one of them on the consensus:

* **Nextclade** with an avian A(H5) dataset, which assigns the clade and
  also reads the cleavage site from the assembled haemagglutinin. Pin
  the dataset tag you used and record it beside the result: the same
  sequence can be placed differently under a later dataset.
* **LABEL**, which IRMA ships, if you already run IRMA. Its H5 module
  classifies the haemagglutinin with hidden Markov models, independently
  of Nextclade. Where the two disagree, the sequence is worth looking at.

The two readings answer different questions and both are worth having.
flucut tells you, while the run is still going and read by read, which
residues the sample carries and in what proportion — including a residue
present in only a part of the population, which a consensus rounds away.
A clade classifier tells you, once a consensus exists, which lineage the
virus belongs to. Neither replaces the other.

Adding a marker
---------------

The table is not in the program sources: it is the metadata file
``metadata/markers.yaml`` of the repository, which the build turns into a
header. A marker is added by writing a record there — its anchor
peptide, where the window lies beside it, the positions with the
convention they are counted in, and one line saying what the residue
means — and rebuilding. The column of ``summary.tsv``, the object of the
JSON summary and the row of the report follow from the record, and no C
source changes.

The technical account, with the anchors, the thresholds and where the
numbers come from, is in ``docs/markers.md`` of the repository.
