The primers
===========

Influenza A and influenza B are amplified by the same multisegment
RT-PCR with different primers, because their segment termini differ.
Which types a run may hold is given by ``--virus`` (``a``, ``b`` or
``both``, the default); it chooses the primers searched and the
references used, and the type of each sample is then read from the
references its reads were located in.

The conserved termini of influenza A
------------------------------------

All influenza A segments end in the same conserved sequences: the 12
nucleotides of one end and the 13 of the other are (almost) identical in
all eight segments and in all subtypes. On the plus strand of an
amplicon (the orientation of the mRNA) they are::

   AGCAAAAGCAGG ... segment ... CCTTGTTTCTACT

The universal primers bind them: Uni12 (``AGCAAAAGCAGG``) at the left
end, Uni13 (``AGTAGAAACAAGG``, the reverse complement of the right end)
at the right end. The fourth base of the 12-nucleotide terminus is A in
some segments and G in others; the degenerate form of Uni12,
``AGCRAAAGCAGG`` (R = A or G), matches both.

The built-in primers of influenza A
-----------------------------------

flucut knows the primers by name; their sequences are quoted from the
publications, each one documented in ``src/seq/primers/primers.h``.

.. list-table::
   :header-rows: 1
   :widths: 18 40 42

   * - name
     - sequence
     - source
   * - uni12
     - ``AGCAAAAGCAGG``
     - the conserved terminus (Uni12 of Hoffmann et al. 2001, as quoted
       by Zhou et al. 2009)
   * - uni12r
     - ``AGCRAAAGCAGG``
     - the same, degenerate fourth base
   * - uni13
     - ``AGTAGAAACAAGG``
     - the conserved terminus (Uni13)
   * - mbtuni-12
     - ``ACGCGTGATCAGCAAAAGCAGG``
     - Zhou et al. 2009, J Virol 83:10309-10313
   * - mbtuni-12.4r
     - ``ACGCGTGATCAGCRAAAGCAGG``
     - Licheri et al. 2025, Genome Med 17:103
   * - mbtuni-13
     - ``ACGCGTGATCAGTAGAAACAAGG``
     - Zhou et al. 2009
   * - uni12-inf1
     - ``GGGGGGAGCAAAAGCAGG``
     - Yang et al. 2025, Microbiol Spectr 13:e02218-24 (there Uni-12/Inf1)
   * - uni12-inf3
     - ``GGGGGGAGCGAAAGCAGG``
     - the same, with G at the fourth base (Uni-12/Inf3)
   * - uni13-inf1
     - ``CGGGTTATTAGTAGAAACAAGG``
     - Yang et al. 2025 (there Uni-13/Inf1)
   * - opti1-f1
     - ``GTTACGCGCCAGCAAAAGCAGG``
     - Mena et al. 2016, eLife 5:e16777 (also written OptiF1)
   * - opti1-f2
     - ``GTTACGCGCCAGCGAAAGCAGG``
     - the same, with G at the fourth base of the terminus (OptiF2)
   * - opti1-r1
     - ``GTTACGCGCCAGTAGAAACAAGG``
     - Mena et al. 2016 (OptiR)

The MBTuni primers carry the 5' tail ``ACGCGTGATCAG``, which overlaps
the terminus by two bases; the Opti1 primers carry ``GTTACGCGCC``,
which does not overlap it. flucut stores each primer whole, with the
position at which its conserved terminus begins, and knows which bases
of a read come from the tail and which from the terminus.

The Opti1 set covers the variable fourth base of the
terminus with two separate primers, where the degenerate ``uni12r`` and
``mbtuni-12.4r`` express it as ``R`` in one.

A protocol which says only "universal primers" may mean the conserved
terminus used bare, as ``uni12`` and ``uni13`` are, or a tailed primer
built on it. Trimming reads of one with the scheme of the other leaves
the tail on them, so the scheme has to match the protocol.

Influenza B: shorter termini, primers of a few segments
-------------------------------------------------------

The segments of influenza B also end in conserved sequences, but
shorter ones: 9 nucleotides at the left end of the plus strand and 10
at the right end. The 10-nucleotide terminus has two forms, which
differ at their sixth base: the HA, the NA and the NS segments have
one, the other five segments the other.

Nine and ten bases are too few for a specific primer, so the universal
primers of influenza B reach five or eight bases further in, into the
untranslated bases next to the terminus. Those bases are not the same
in every segment: they are shared by a few segments at a time, which is
why the primer set of influenza B is a cocktail of six pairs for the
eight segments. Each primer carries, in flucut, the segments it was
designed for:

====================== ====================
primer pair            segments
====================== ====================
b-pbs-unif/-unir       PB1 and PB2
b-pa-unif/-unir        PA
b-hana-unif/-unir      HA and NA
b-np-unif/-unir        NP
b-m-uni3f, b-mg-uni3f, M
b-m-uni3r
b-ns-uni3f/-uni3r      NS
====================== ====================

They are the universal IBV-GA2 primer cocktail of Zhou et al. 2014,
J Clin Microbiol 52:1330-1337 (Table 1); the two forms of the primer of
the left end of the M segment are both in the cocktail. Each primer is
a 5' tail (``GGGGGG`` at the left end, ``CCGGGTTATT`` at the right end,
neither of them influenza sequence) followed by the terminus of its end
and by the bases next to it.

This is why what flucut keeps at a segment end depends on the primer
found there. With references, the bases of a primer are kept at an end
only when the read was located in a segment that primer binds: there
they are the sequence of the segment. At an end of any other segment
they are the primer's own sequence, as at an internally primed end, and
the whole primer is removed. The universal primers of influenza A bind
every segment, so for them this never applies.

Keeping only the 9 or 10 universal bases at such an end is not
possible: the bases which differ between the segments lie between the
template and the universal terminus, so what is kept at an end is
either the whole stretch of segment sequence under the primer or none
of it.

Influenza B has no subtypes. Its diversity is two lineages, Victoria
and Yamagata, which the records of the sequence databases rarely name;
flucut therefore reports a sample of influenza B as ``B``, and as
``B/Victoria`` or ``B/Yamagata`` only where the description of the
reference chosen says which lineage it is. The lineage of a sample is
better called downstream, from its assembled HA.

.. _influenza-b-features:

What flucut does for influenza B
--------------------------------

The trimming and everything built on it work the same for both types:
the primers are found at both read ends, the bases which are not the
sequence of the segment are removed (with the rule above at the ends of
influenza B), the read ends are located on the references and the reads
are put in their categories (full template, internal priming, internal
deletion, fragment, unplaced, unknown), written by category with
``--by-category``, and counted in the reports, with the segments, their
coverage and the internal priming sites. The HA cleavage site and the
fusion peptide are read as for influenza A: the anchor ``G..GAIAG``
matches the ``GFFGAIAG`` of influenza B (:doc:`cleavage_site`).

These are of influenza A only:

- **the subtype**: a sample of influenza B is ``B``, or ``B/Victoria``
  and ``B/Yamagata`` where the reference chosen names the lineage
  (above);
- **the reads of a second HA or NA subtype**: not counted for influenza
  B, whose lineages are not subtypes (:doc:`reports`);
- **the markers** (:doc:`markers`): no anchor of the table is found in
  influenza B, so no marker is reported;
- **the choice of the scheme from the tails** (``-s auto``, below):
  it chooses among the schemes of influenza A, and the primers of
  influenza B are added to the scheme chosen;
- **flucut-posthoc**, the check of the read categories against
  mappings (``posthoc/``), searches the termini of influenza A and
  numbers the segments as influenza A does.

Influenza B has been tested on reads simulated from three genomes (of
the Victoria and the Yamagata lineage, and one of 1954, before the two
diverged) and on four public Illumina runs of B/Victoria, paired reads of
fragmented amplicons of the IBV-GA2 primers (BioProject PRJNA1426172).
No Nanopore run of influenza B has been tested.

The schemes
-----------

A scheme names the primers expected at each end of an amplicon
(option ``-s``):

======================= ============================ ====================
scheme                  left end                     right end
======================= ============================ ====================
uni12-mbt13+ibv-ga2     uni12r and the seven of      mbtuni-13 and the
                        influenza B                  six of influenza B
uni12-mbt13             uni12r                       mbtuni-13
ibv-ga2                 the seven of influenza B     the six of
                                                     influenza B
mbt                     mbtuni-12, mbtuni-12.4r      mbtuni-13
uni                     uni12, uni12r                uni13
opti                    opti1-f1, opti1-f2           opti1-r1
uni-inf                 uni12-inf1, uni12-inf3       uni13-inf1
======================= ============================ ====================

The default, ``-s auto``, chooses the scheme of influenza A from the
reads (see below), and adds the primers of influenza B when ``--virus``
includes it. Where the reads do not single out a scheme, the scheme of
``--virus`` is used: ``uni12-mbt13+ibv-ga2`` for both types,
``uni12-mbt13`` for influenza A alone (the untailed degenerate Uni12
primer and the tailed Uni13 primer, the combination used at the
Institut Pasteur du Cambodge); with ``--virus b`` it is always
``ibv-ga2``, the only scheme of influenza B. When a scheme has several primers at an end, the primer whose
match reaches furthest out is taken, then the one with fewest errors.
Several schemes joined by ``+`` are one scheme with the primers of each,
e.g. ``-s opti+ibv-ga2`` for a run of the Opti1 primers which may also
hold influenza B.

A scheme of both types costs almost nothing: once the references of a
sample say which type it is, only the primers of that type are searched
in its reads, and the trimmed reads of an influenza A sample are the
same as they were before influenza B was supported.

What it does cost is fixed, once per run: the references of influenza B
in the default set, and the automatic error rate (``-e auto``), which is
measured on the first 1000 reads before the type is known and therefore
searches every primer of the scheme. Together they were 0.9 s on a run
of 183000 Nanopore reads of influenza A (6.11 s against 5.25 s before
influenza B was supported, with the default references of each). Being
fixed, they are a smaller and smaller share of a longer run: on 3.2 GB
of reads in four barcodes, 182.2 s against 181.4 s. ``--virus a`` or
``--virus b`` on a run which holds one type avoids them.

The rate itself is measured per type: at the few read ends where a
primer of the other type matches best, it matches with many errors,
and counting those with the rest would raise the rate and with it the
errors allowed in a primer (on that run, 4.33% instead of 3.69%, which
is on the other side of the threshold at which 30% is allowed rather
than 20%). The rate reported is that of the type whose primers cover
most of the read ends.

Choosing the scheme from the reads
----------------------------------

``-s auto``, the default, chooses the scheme of influenza A from the
first reads of the first input file. The schemes differ in the tail in front of the
conserved terminus, so the tails are what is counted:

#. At the start of each read and of its reverse complement (within the
   window of the primers, ``-w``), the conserved terminus of each end
   is searched, with at most one error. The termini of influenza B
   differ from those of influenza A at two or three bases, so they are
   not counted.
#. Where a terminus is found, each tail of the built-in schemes of
   influenza A for that end (``ACGCGTGATC`` of ``mbt``, ``GTTACGCGCC``
   of ``opti``, ``GGGGGG`` and ``CGGGTTATT`` of ``uni-inf``) is aligned
   to the bases directly in front of it, with at most 20% of its length
   in errors. A tail elsewhere near the end of a read does not count: a
   run of G at the end of a read of two-colour Illumina chemistry, where
   no signal is read as G, is not a primer.
#. Each end is then judged by the share of the read ends with the
   terminus which carry a tail:

   - **tailed**, if the most frequent tail is in at least 5% of them
     and in at least ten times as many as the next tail;
   - **two tails**, if it has the share but not the lead over the next
     one: a run holding more than one protocol;
   - **untailed**, if no tail has 5%;
   - **too few**, if fewer than 100 read ends have the terminus.

   The share is over the read ends with the terminus, not over all the
   reads: on Illumina reads only a fragment starting at an amplicon end
   has a terminus.
#. The scheme chosen is the one whose tails are those of the two ends
   (an end with too few read ends fits any scheme), if at least one end
   has a tail.

The first 2000 reads are counted; unless both ends are then tailed, the
reads up to the first 20000 are counted before deciding. On the public
runs of the three protocols the tail of the protocol was in front of
81% to 98% of the read ends with a terminus at each end, and in runs
without a tail no tail was in front of more than 0.6% of them. The 5%
share is low on purpose: some runs mix a tailed and an untailed primer
at an end, and the share of the tail is then anywhere between 13% and
66%. A higher threshold would split such runs between two schemes; with
5% they are all seen as having the tail.

Where no scheme is chosen, the default scheme of influenza A
(``uni12-mbt13``) is used, and the reason is reported:

- **no tail**: neither end has one. The scheme without tails (``uni``)
  would match the reads, but it is not chosen: its bare termini, twelve
  and thirteen bases, match inside the template by chance more often
  than the tailed primer of the default, and on simulated reads of
  untailed primers it removed template bases which the default kept;
- **two tails at an end**: the run holds more than one protocol, and
  any single scheme would be wrong for part of it;
- **too few read ends with a conserved terminus**: on libraries of this
  family of protocols this usually means that the reads are not
  amplicons of the universal primers at all (for example random cDNA);
- **no built-in scheme has the tails found**: for example an untailed
  left end and the tail ``CGGGTTATT`` at the right end.

The choice is within influenza A. With ``--virus both`` the primers of
influenza B (``ibv-ga2``) are added to the scheme chosen, and the primers
of the type of each sample are searched in its reads once its references
say which type it is, as with any scheme of both types. With ``--virus
b`` there is nothing to choose.

The scheme chosen, the outcome and the counts of each tail at each end
are in the messages, in the settings of the HTML report and in the JSON
summary (``scheme.choice``). In the live mode, the scheme is chosen from
the first file ready, and the reports are written from then on.

The choice is a fixed cost, once per run, beside those of a scheme of
both types (above): 0.08 s of processor time on a run of the Opti1 primers, where
both ends carry a tail and 2000 reads decide, and 0.6 s on the run of
171086 Nanopore reads of the default scheme, whose untailed left end
makes it read 20000 (the median of seven repetitions, on a shared
machine under load; ``-D`` prints it).
