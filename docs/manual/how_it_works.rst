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

The whole primer and its shorter forms are searched in one alignment,
in which the missing outer bases cost as many errors as the shorter
form is allowed fewer: the match with the fewest errors, counting this
cost, is taken, and among equally good matches the outermost. The
overhang is considered only if the primer is not found inside the read.

The window
----------

A match of a primer deep inside a read is not the primer which started
the amplicon: trimming there would remove template. The primer of an
end is therefore accepted only if it begins (left end) or ends (right
end) within ``-w`` bases of that end of the read (default 100: enough
for the barcoding adapters of Oxford Nanopore, 46 bases, and further
flanking bases).

Within the window, the outermost of the best matches of a primer is
taken, and the outermost of the matches of the primers of an end (e.g.
the two forms of a degenerate primer): the primer which
started the amplicon is the one nearest to the end of the read, and a
primer-like sequence further inside is template. Matches of a primer
between the two end primers are counted (tag ``ni``), not trimmed: they
are the sites at which a primer can anneal, and the junctions of
chimeras.

The adapters
------------

At an end of the read at which no primer is found (the molecule was
broken, or the read begins too far inside the primer), flucut searches
the sequencing adapters of Oxford Nanopore (``--adapters``, default: all
the kits): the native barcoding adapters (``nb``, e.g. SQK-NBD114), the
rapid barcoding adapter (``rbk``, e.g. SQK-RBK114) and the ligation
adapter (``lsk``, e.g. SQK-LSK114), with the sequences of Dorado. The
barcode is searched as 24 Ns, which match any base: any barcode of the
kit is recognised, none is identified (the reads normally come already
demultiplexed). The errors allowed are 20% of the fixed bases of an
adapter (e.g. 4 of the 22 of a native barcoding adapter); an adapter
without a barcode may overhang the read by up to half of its length,
one with a barcode may not (the part left in the read could then be
mostly the barcode, which matches anything). The adapter must lie in
the window of its end; of the adapters found, the innermost is taken
(the barcode adapter, not the ligation adapter before it).

An adapter found is removed with everything outside it, and the part of
the read between the adapters is analyzed again: a primer of which only
a few bases are left is now at an end, next to the adapter, where it can
overhang (up to half of it missing), and is removed as well. In this
second analysis, the primer must lie within 5 bases of the adapter: in
the usual window of the new end, a primer-like sequence of the template
could otherwise be taken for it. The tags ``a5`` and ``a3`` tell which
kit's adapter was removed at the start and at the end of a read.

The search uses the adaptor search of Lumipore (``src/adaptor-search.c``,
unchanged); the adapters are searched only at the ends without a primer,
so the reads with both primers are trimmed as without it.

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

The files are trimmed by ``-j`` worker threads, started for each file
(as in Lumipore): each thread reads a batch of 1000 reads in turn,
trims it and compresses it with its own objects, and writes it in the
order of the input, waiting for its turn. Each batch is a gzip member of
its own, so that the output file is the same, byte for byte, for any
number of threads. Each thread counts its reads; the counts are summed
at the end of the file.
