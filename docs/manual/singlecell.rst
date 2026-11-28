.. _singlecell:

Single cell
-----------

The ``singlecell`` directory contains barcode whitelists and correction,
cell calling, and library complexity modules. These tools support a
common two-pass workflow: identify cells from exact barcode counts, then
correct barcodes using only the selected cells. For edit distance over
short candidate lists, see :ref:`barcodes`.

.. _whitelist:

Barcode whitelists (``whitelist.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A whitelist is the set of the possible barcodes of a tag: e.g. the cell
barcodes of a single-cell platform (millions), the barcodes of the cells
of an experiment (thousands), or the barcodes of a small multiplexing
kit. The barcodes have all the same length, up to 32 bases, and are
stored encoded (2 bits per base, in a ``uint64_t``) in a hash set.

**Reading a whitelist.** ``whitelist_new(path, logger)`` reads a text
file (optionally gzipped) with one barcode per line: only the first word
of each line is used, without a suffix starting with ``-`` (such as the
``-1`` of the barcode files of Cell Ranger), so that e.g. a
``barcodes.tsv.gz`` file of Cell Ranger can be used directly. Empty lines
are skipped; the barcodes must consist of the letters ``ACGT``. Invalid
files exit the program with a message giving the line::

   ERROR | Whitelist cells.txt, line 3: barcodes of different lengths

``whitelist_new_from_codes(codes, n, length)`` creates a whitelist from
encoded barcodes (e.g. the cells selected by the cell caller).

**Correcting a barcode.** ``whitelist_match()`` searches a barcode at one
end of a sequence, the *anchored* end (e.g. the end next to an adaptor):
the barcode of the whitelist is aligned completely, while the sequence
can be aligned partially, as the barcode is followed by other elements
(e.g. a UMI). The barcode with the smallest edit distance (substitutions,
insertions and deletions), up to a maximum (at most
``WHITELIST_MAX_DIST``, 3), is returned if it is unique:

.. code-block:: c

   whitelist_match_t m;
   /* the barcode at the start of seq (from_end false), up to 2 edits */
   whitelist_match(wl, seq, len, false, 2, &m);
   switch (m.status) {
     case WHITELIST_EXACT:       /* found with distance 0 */
     case WHITELIST_CORRECTED:   /* unique closest barcode, m.dist > 0 */
       whitelist_decode(wl, m.barcode, bc);   /* length + 1 bytes */
       umi = seq + m.consumed;   /* the bases after the barcode */
       break;
     case WHITELIST_AMBIGUOUS:   /* several barcodes at the minimum
                                    distance */
     case WHITELIST_NOT_FOUND:   /* none within the maximum distance */
       break;
   }

``m.consumed`` is the length of the aligned part of the sequence: with an
insertion or a deletion in the barcode, it differs from the barcode
length, and the following elements start at a different position. With
``from_end`` true, the barcode is searched at the end of the sequence
(e.g. before an adaptor).

The search enumerates, for each distance, all the sequences at that
distance from the read and looks them up in the set: the time thus grows
quickly with the distance, and a maximum distance of 2 is the practical
limit for barcodes of 16 bases. With a large whitelist (e.g. all the
barcodes of a platform), a barcode with errors is often closer to another
barcode of the whitelist than to the correct one: correct the barcodes
using only the barcodes of the cells (see below).

Other functions: ``whitelist_encode()`` and ``whitelist_decode()``
(between sequences and codes), ``whitelist_contains()``,
``whitelist_size()``, ``whitelist_barcode_length()`` and
``whitelist_codes()`` (all the codes, in an allocated array). A
whitelist can be shared by several threads after its construction.

.. _cell-caller:

Identifying the cells (``cell-caller.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In a single-cell library, most barcodes of the whitelist of the platform
are not used, and the reads of the cells are mixed with background reads
(e.g. of empty droplets). The cell caller identifies the cells from the
counts of the barcodes found **exactly** in the reads, given the
expected number of cells *N*, with the rule of Cell Ranger 2: the cells
are the barcodes with a count of at least 1/10 of the 99th percentile of
the counts of the *N* most frequent barcodes (the threshold is rounded to
the nearest integer, as in Cell Ranger, and is at least 1).

.. code-block:: c

   cell_caller_t *cc = cell_caller_new();
   /* first pass over the reads: count the exact barcodes */
   whitelist_match(platform_wl, seq, len, false, 0, &m);
   if (m.status == WHITELIST_EXACT)
     cell_caller_add(cc, m.barcode);
   ...
   whitelist_t *cells = cell_caller_run(cc, 5000,
       whitelist_barcode_length(platform_wl), NULL,
       "out/cell_barcodes.tsv", logger);
   cell_caller_delete(cc);
   /* second pass: correct the barcodes using only the cells */
   whitelist_match(cells, seq, len, false, 2, &m);

``cell_caller_run()`` returns the cells as a new whitelist; if a path is
given, the cells are written to it, with their counts (``barcode<TAB>
count``). Cells known from another source (e.g. the cells called by Cell
Ranger in an Illumina run of the same library) can be given as a
whitelist (``known``): they are added to the selected cells (useful for
targeted libraries, with few reads per cell). In a multithreaded
program, each thread counts with its own cell caller, and the callers are
merged (``cell_caller_merge()``) before the selection.
``cell_caller_sorted_counts()`` returns the counts in decreasing order,
e.g. for plotting the "knee plot" of the barcodes.

.. _complexity:

Library complexity (``complexity.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The complexity of a library is the number of its distinct molecules. It
tells whether sequencing it more deeply is worthwhile: if most reads are
duplicates, more reads will mostly be more duplicates. Without an
alignment, the copies of a molecule are recognized by a *key*, e.g. the
UMI and the first and last bases of the insert; the distinct keys are
counted, e.g. with a HyperLogLog sketch (:ref:`ds`).

Sequencing errors make the copies of a molecule look distinct. With an
error rate *e* per base (e.g. measured on the alignments of known
adaptor sequences), a key of length *L* is without errors with
probability *p* = (1 − *e*)\ :sup:`L`; each read with an error in the key
is assumed to add a distinct key. The *pN* reads with correct keys thus
contain *D* = *D*\ :sub:`obs` − (1 − *p*)\ *N* distinct molecules, and the
number of molecules *M* of the library is estimated, as in Picard
``EstimateLibraryComplexity``, by solving the Lander-Waterman equation
*D* = *M* (1 − exp(−*pN*/*M*)).

.. code-block:: c

   complexity_t c;
   complexity_estimate(n_reads, hll_estimate(&keys), mean_key_length,
                       error_rate, &c);
   if (!c.estimated)
     ;   /* fewer than COMPLEXITY_MIN_READS reads with a correct key */
   else if (c.saturated_unknown)
     ;   /* no duplicates: only a lower bound (c.molecules_seen) */
   else
     printf("duplicates: %.1f%%, molecules: %.0f, with 2x reads: %.0f\n",
            100 * c.duplication, c.molecules, c.molecules_2x);

The fields of the result: ``error_free`` (*p*), ``molecules_seen`` (the
distinct molecules in the reads), ``duplication`` (the fraction of
duplicate reads), ``molecules`` (*M*) and ``molecules_2x`` (the distinct
molecules expected with twice the reads). The estimate assumes molecules
sampled uniformly (no amplification bias) and independent errors; with
few duplicates, it is sensitive to the error of the count of the
distinct keys. In simulations with 1–3% errors, the
duplication was within about 0.07 of the true value, and the number of
molecules within about 30% for libraries with more than 20% duplicates.
