.. _barcodes:

Barcode and short-sequence distances
------------------------------------

The sequence component provides a bounded edit distance for comparing
short sequences with a small list of candidates. Single-cell barcode
whitelists, cell calling and library complexity are documented in
:ref:`singlecell`.

.. _edit-distance:

Edit distance of short sequences (``edit-distance.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The single-cell :ref:`whitelist` module finds a matching barcode by
enumerating the sequences around the read, which is fast for a large set
of barcodes and a small distance. When a sequence is instead compared
with each sequence of a short list (e.g. the indices of the samples of a
sequencing run, which can have different lengths), the edit distance is
computed directly.
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
