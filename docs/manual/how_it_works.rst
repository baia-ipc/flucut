How it works
============

.. figure:: figures/read_layout.*
   :width: 100%

   The parts of a read of an amplicon, and what each trimming policy
   keeps (``make_diagrams.py``).

The layout of a read
--------------------

A read of an amplicon contains, from one end to the other: flanking
bases (adapter, barcode, or bases of the library preparation), the
primer of one end, the template (the segment), the primer of the other
end, and flanking bases again. It can be sequenced in either
orientation, and it can begin or end inside a primer.

flucut searches the primers of the left end of the scheme at the start
of the read and those of the right end at its end. If the read was
sequenced in the other orientation, the primers are found in its reverse
complement; the positions are always reported in the original read,
which is not reverse-complemented.

The search of a primer
----------------------

A primer is aligned to the read with edit distance (substitutions,
insertions and deletions), allowing ``-e`` errors (default 0.2 of the
length of the primer: 2 errors in Uni12). Degenerate bases (e.g. the R
of Uni12) match any of the bases they stand for, in the primer and in
the read (e.g. N).

At an end of the read, the primer can be partial:

- it can overhang the read, when the read begins inside the primer: the
  missing part is not counted as errors, and up to half of the primer
  can be missing;
- it can lack up to half of its outer bases after some flanking bases,
  when the amplicon was shortened before the adapter was ligated; this
  shorter form must match with one error less than its length allows,
  as a shorter sequence matches by chance more easily.

The window
----------

A match of a primer deep inside a read is not the primer which started
the amplicon: trimming there would remove template. The primer of an
end is therefore accepted only if it begins (left end) or ends (right
end) within ``-w`` bases of that end of the read (default 100: enough
for the barcoding adapters of Oxford Nanopore, 46 bases, and further
flanking bases).

Within the window, the outermost match is taken: the primer which
started the amplicon is the one nearest to the end of the read, and a
primer-like sequence further inside is template. Matches of a primer
between the two end primers are counted (tag ``ni``), not trimmed: they
are the sites at which a primer can anneal, and the junctions of
chimeras.

The trimming
------------

The removal is by position: the whole span of the primer alignment is
removed, whether or not the bases under it differ from the primer, and
everything outside it.

``--trim primer`` (the default)
   removes every primer-derived base, tail and conserved terminus. No
   base of the output comes from a primer. The consensus is shorter than
   the segment by the termini (12 and 13 bases): those bases would have
   been the primer's, not measured sequence.

``--trim terminus``
   removes the tail and everything outside the primer, and keeps the
   bases of the conserved terminus, so that the consensus has complete
   ends. The kept bases are those of the primer: with the degenerate
   Uni12, the fourth base is a mix of A and G; where the segment differs
   from the primer (e.g. the 3' end of segment 2 of A/Puerto Rico/8/1934,
   ``CCTTGTTCCTACT``), the primer's base is kept. At an internal priming
   site, the primer bases are not the template at all: this policy
   therefore requires the references (``-r``), against which such ends
   are recognised and trimmed whole (:doc:`internal_priming`).

Reads with nothing left (primer dimers) and reads shorter than ``-m``
bases after the trimming are not written.

Threads
-------

The reads are read in batches; each batch is split between the threads
(``-j``), each of which trims and compresses its part, and the parts are
written in the order of the input: the output is the same for any
number of threads.
