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
publications, each one documented in ``src/seq/primers.h``.

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

The tailed primers carry the 5' tail ``ACGCGTGATCAG``, which overlaps
the terminus by two bases: flucut stores each primer whole, with the
position at which its conserved terminus begins, and knows which bases
of a read come from the tail and which from the terminus.

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
======================= ============================ ====================

The default is the scheme of ``--virus``: ``uni12-mbt13+ibv-ga2`` for
both types, ``uni12-mbt13`` for influenza A alone (the untailed
degenerate Uni12 primer and the tailed Uni13 primer, the combination
used at the Institut Pasteur du Cambodge) and ``ibv-ga2`` for influenza
B alone. When a scheme has several primers at an end, the primer whose
match reaches furthest out is taken, then the one with fewest errors.

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
