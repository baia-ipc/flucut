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

Reads with several amplicons
----------------------------

Some reads hold two or more amplicons, joined end to end: molecules
ligated to each other in the library preparation, or sequenced one after
the other without the read being split. The primer of an inner amplicon
is then inside the read, not at its end. If that amplicon was internally
primed, its primer bases would stay in the read and reach the consensus:
in a sample with few complete HA amplicons, such reads alone can be
enough to give AAA at the site.

flucut therefore splits these reads at their junctions, and trims each
part as a read of its own (:doc:`output`). A junction is accepted in two
cases:

- the end of an amplicon followed, within 100 bases, by the start of
  another: e.g. the reverse complement of Uni13, the adapters and
  barcodes, then Uni12. This needs no reference;
- with references, a single primer inside the read at which the read
  stops following the segment: the bases on one side of it are located
  in a reference, and those on the other side are not found next to
  them.

A read of a complete amplicon which runs through a primer-like site of
its segment is not cut: the segment continues on both sides of the site.
A primer of the same end found at both ends of a read (two amplicons
joined without a primer between them) is removed at both ends.

``summary.tsv`` reports the reads split (``split``); ``--no-split`` turns
the splitting off. Without references (``-r off``) only the first kind
of junction is found.

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

With references, an end classified as internal (or whose origin cannot be
determined: cut when unsure) loses its whole primer, and only an end
confidently classified as terminal keeps its conserved terminus. With
``-r off`` the classification is not done and every primer is removed.

The built-in references
-----------------------

By default (no ``-r``) flucut classifies the reads against a **built-in
reference set**: complete segments, with the UTRs, of several influenza A
subtypes (H1N1, H3N2, H5N1, H7N9, H9N2). It is not shipped with the
program; it is downloaded from NCBI on the first run that needs it and
cached (under ``$XDG_CACHE_HOME/flucut`` or ``~/.cache/flucut``), so later
runs do not download it again. If it cannot be downloaded (no network, no
``curl``), flucut warns and runs without references.

Each barcode (a different sample) is classified against one reference
of each segment: the one closest to its reads. Each segment is chosen
by its own reads, so that a reassortant gets, for each segment, the
reference of the subtype that segment is closest to; a segment whose
reads do not single out a reference follows the clear ones. This matters
because the UTRs are conserved across subtypes (so a terminal read
matches the segment of any subtype) and the coding regions diverge (so
a read with a noisy start finds spurious interior matches); using one
reference per segment avoids scattering the reads and inflating the
internal priming.

The choice is made from the first 2000 reads of the first file of the
barcode. Eight probes are taken along each read and located in the
references. A probe votes for the reference it is located in, unless
another reference matches it as well: such a probe tells nothing, and
does not vote.

A segment is clear if its reference with most votes has at least 20
votes and at least 60% of the votes of the segment; that reference is
chosen. The other segments follow the clear ones where they can: a
segment with too few votes takes the reference of the subtype chosen in
most clear segments, and so does a segment without a clear lead, if
that reference has at least half the votes of the leading one.
Otherwise a segment without a clear lead takes its reference with most
votes (common for the internal genes, which differ little between
subtypes), and for a segment with too few votes the reads up to 20,000
are used as well; if it still has too few, it keeps all its references.
The messages (``-v``) give the reference chosen for each segment, with
its votes and how it was chosen.

For a run of one segment alone (e.g. an HA amplicon), ``--segment``
does the same for the segments: the read ends are classified against
the references of that segment only (``--segment HA``), which is faster
and places no read end in another segment by a chance match. On reads
of the whole genome it must not be used, as the ends of the other
segments are then not located (:doc:`cli_reference`).

``-r <file>`` uses your own references instead of the built-in set (e.g.
the exact segments of the run's isolate), and ``-r off`` turns references
off. ``scripts/fetch_references.py`` downloads segment references from
NCBI, with their provenance.

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
