The primers
===========

The conserved termini
---------------------

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

The built-in primers
--------------------

flucut knows the primers by name; their sequences are quoted from the
publications, each one documented in ``src/primers.h``.

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

The schemes
-----------

A scheme names the primers expected at each end of an amplicon
(option ``-s``):

=============== ============================ ====================
scheme          left end (Uni12 side)        right end (Uni13)
=============== ============================ ====================
uni12-mbt13     uni12r                       mbtuni-13
mbt             mbtuni-12, mbtuni-12.4r      mbtuni-13
uni             uni12, uni12r                uni13
=============== ============================ ====================

``uni12-mbt13`` is the default: the untailed degenerate Uni12 primer and
the tailed Uni13 primer, the combination used at the Institut Pasteur du
Cambodge. When a scheme has several primers at an end, the best match is
taken.
