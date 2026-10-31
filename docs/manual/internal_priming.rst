Internal priming
================

The artefact
------------

.. figure:: figures/internal_priming.*
   :width: 100%

   The Uni12-like site of the HA of A/chicken/Cambodia/X0815302/2013
   (H5N1), NCBI KF918483.1, and the amplicon the primer starts there
   (``make_diagrams.py``).

The untailed Uni12 primer is twelve bases long, with a degenerate base.
A sequence similar enough to it can occur inside a segment: the primer
anneals there and starts an amplicon in the middle of the segment. The
reads of such an amplicon begin with the bases of the primer where the
template has others, so at that position the reads carry the primer's
sequence, not the segment's. If enough reads are primed at the same
site, the consensus shows the primer's base there: a variant which looks
real and is not. Such sites recur from run to run in a given set of
sequences, so the false variant is reproducible, which makes it easy to
mistake for a true one.

What flucut does about it
-------------------------

The primer bases of an internally primed read are at its end, like those
of any other read: removing every primer (the default, without
references) removes the artefact, without any reference. The site is still
covered by the reads of the complete amplicons, which carry the
template.

In the benchmark (:doc:`performance`), the untrimmed simulated reads give
AAA at the site in the IRMA consensus, and the reads trimmed by flucut
give ATA, the truth.

Classification against references
---------------------------------

With reference sequences of the segments (``-r``, a FASTA file), flucut
also tells where each read started. The bases which follow the primer of
an end (the probe: ``--probe-len``, default 32 bases; for the right end,
the reverse complement of the bases before its primer) are located in
the references:

- at a terminus of a segment (within ``--terminal-window`` bases of it,
  default 60, as references often lack some bases of the termini): a
  normal amplicon (tags ``ol``/``or``: ``terminal``);
- at an interior position: internal priming (``internal``); the tags
  ``ml``/``mr`` give the reference and the position;
- not located (no reference of the subtype, a probe too short):
  ``unknown``.

Among equally good locations, the one at a terminus is taken: an
ambiguous read is never classified as internally primed.

The references are indexed by their k-mers (10-mers), so that the time
hardly depends on their number; a probe without a candidate is aligned
to all of them. ``--exhaustive-search`` aligns every probe to all of
them, to check the index (the outputs were identical in the benchmark).

Without references the classification is not done and every primer is
removed. With references, an end classified as internal (or whose origin
cannot be determined: cut when unsure) loses its whole primer, and only
an end confidently classified as terminal keeps its conserved
terminus.

``scripts/fetch_references.py`` downloads segment references from NCBI,
with their provenance.

The references must be **complete segments, with the untranslated
regions (UTRs)**, not coding sequences (CDS). A read is located by the
sequence just after its 5' primer, which in a full segment is the 5'
UTR; a CDS reference (starting at the ATG, ending at the stop codon)
does not contain it, so the full-length reads cannot be located and only
the internally primed reads, whose probe falls inside the ORF, map —
which would make the classification and the coverage plot misleading.
flucut warns when the references look like coding sequences
(:doc:`troubleshooting`); a complete segment begins with the conserved
terminus (``AGCRAAAGCAGG``), not with ``ATG``.

The length check
----------------

The complete amplicons of a segment (both primers found) have the length
of the segment. Internal priming, or a deletion (e.g. in the defective
interfering RNAs of the polymerase segments), shows as a group of
shorter amplicons. flucut writes their lengths to ``lengths.tsv``
(:doc:`output`):

- with references, for each reference the amplicons shorter than 90% of
  its length are counted (``short_complete`` in ``summary.tsv``) and
  reported, with the length where most of them are. The reference
  length is the measure, not the most frequent length, as the
  internally primed amplicons can be the majority;
- without references, the peaks of the lengths are listed.

In the benchmark, 209 of the 340 complete HA amplicons were reported as
short, most at 670-679 bases: the length from the internal site to the
end of the segment.
