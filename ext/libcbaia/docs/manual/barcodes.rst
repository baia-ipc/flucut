.. _barcodes:

Barcodes, cells and library complexity
--------------------------------------

Four modules of the ``seq`` component implement analyses which are
common to the tools working with barcoded reads: the correction of
barcodes using a whitelist, the edit distance of short sequences (e.g.
for matching indices to a list), the identification of the cells of a
single-cell library, and the estimate of the complexity of a library.

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

.. _edit-distance:

Edit distance of short sequences (``edit-distance.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The whitelist finds the closest barcode by enumerating the sequences
around the read, which is fast for a large set of barcodes and a small
distance. When a sequence is instead compared with each sequence of a
short list (e.g. the indices of the samples of a sequencing run, which
can have different lengths), the edit distance is computed directly.
``edit-distance.h`` computes it for short sequences (up to
``SHORT_EDIT_DISTANCE_MAX``, 64 bases), **bounded** by a maximum
distance *k*: the dynamic programming is limited to the diagonal band of
width 2\ *k* + 1, and stops as soon as all the values of a row exceed
*k*, so that the comparison of dissimilar sequences ends after a few
rows. The functions are ``static inline``: the comparison with all the
sequences of a list is often in the innermost loop of a program, and
is compiled in place, without a function call. For longer sequences, an
alignment library (e.g. edlib) is used by the programs.

``long short_edit_distance(const char *a, size_t alen, const char *b, size_t blen, long k)``
   the edit distance (substitutions, insertions and deletions) of ``a``
   and ``b``, if at most ``k`` (``k`` >= 0), otherwise ``k + 1``;
``long short_prefix_distance(const char *query, size_t qlen, const char *target, size_t tlen, long k, size_t *end)``
   the distance of the best alignment of the whole ``query`` with a
   prefix of ``target`` (not empty), if at most ``k``, otherwise -1;
   ``*end`` is set to the length of that prefix (the shortest among those
   with the best distance, as edlib in the mode ``SHW``). This is the
   search of an element at the start of a read, followed by other
   elements: e.g. ``short_prefix_distance("ACGT", 4, "AGTTTT", 6, 1,
   &end)`` returns 1 (the deletion of ``C``), with ``end`` 3.

For example, the closest index of a list, if unique:

.. code-block:: c

   /* the index of the read matched to the indices of the samples:
      the sample with the smallest distance, up to max_dist; -1 if none
      or ambiguous */
   long best = max_dist + 1, sample = -1;
   for (size_t s = 0; s < n_samples; s++) {
     long d = short_edit_distance(index, index_len, samples[s].index,
                                  samples[s].index_len, best);
     if (d < best) {
       best = d;
       sample = (long)s;
     } else if (d == best && d <= max_dist)
       sample = -1;   /* another sample at the same distance */
   }

(passing the best distance found so far as the bound makes the
comparisons with the following indices faster).

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
distinct keys. In the simulations of Lumipore, with 1–3% errors, the
duplication was within about 0.07 of the true value, and the number of
molecules within about 30% for libraries with more than 20% duplicates.
