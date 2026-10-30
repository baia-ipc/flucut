.. _ds:

Data structures
---------------

The ``ds`` component contains data structures, which only depend on the
``core`` component. Hash tables are provided by khash (``khash.h``, in
``ext/khash``), included with the library and available to the projects
(see :ref:`khash`).

Bit operations (``bitwise.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Macros for single bits of an integer of any width: ``GET_BIT(e, n)``,
``SET_BIT(e, n)``, ``UNSET_BIT(e, n)`` (the bit number can be up to the
width of the type of ``e`` minus 1).

``bitwise_mix64(x)`` mixes the bits of a 64-bit integer, so that each bit
of the result depends on all bits of ``x`` (the finalizer of splitmix64):
e.g. to hash an integer key, as the Bloom filter does, or to improve a
weak hash value, as ``hll_hash()`` does.

Bloom filter (``bloom-filter.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A Bloom filter of integer keys (``size_t``), of 2\ :sup:`k` bytes, with
two hash functions derived from the splitmix64 hash of the key. A query
returns ``false`` if the key was certainly not inserted, ``true`` if it
may have been:

.. code-block:: c

   bloom_filter_t *bf = bloom_filter_new(20);   /* 1 MiB */
   bloom_filter_insert(bf, key);
   if (bloom_filter_query(bf, other_key))
     ...   /* possibly present: check in the exact structure */
   bloom_filter_delete(bf);

Set of integers (``set.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A set of ``set_elem_t`` (``size_t``) values, for large numbers of
elements with little memory overhead: a sorted array, queried by binary
search, plus an unsorted insertion buffer, which is sorted and merged
into the array when full. Optionally, the elements are divided into
2\ :sup:`k` buckets (by their lowest *k* bits), and a Bloom filter speeds
up the queries of absent elements.

.. code-block:: c

   /* buffer of 4096 elements, Bloom filter of 2^16 bytes, no buckets */
   Set *set = set_new(4096, 16, 0);
   for (...)
     set_insert(set, value);        /* set_insert_simple(): no check of
                                       duplicates */
   set_finalize_insertions(set);    /* faster queries afterwards */
   if (set_query(set, x))
     ...
   set_delete(set);

Heap (``heap.h``)
~~~~~~~~~~~~~~~~~

A binary min-heap of pointers, ordered by a comparison function with the
semantics of ``qsort()`` (called with the stored pointers): a priority
queue, or the merging of sorted sequences (as in the file merger).

.. code-block:: c

   static int cmp_ints(const void *a, const void *b) {
     int x = *(const int *)a, y = *(const int *)b;
     return (x > y) - (x < y);
   }

   heap_t *h = heap_new(cmp_ints, 16);   /* initial capacity; it grows */
   int values[] = {5, 1, 3};
   for (size_t i = 0; i < 3; i++)
     heap_push(h, &values[i]);
   while (!heap_is_empty(h))
     printf("%d\n", *(int *)heap_pop(h));   /* 1, 3, 5 */
   heap_delete(h);                           /* not the elements */

HyperLogLog (``hyperloglog.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The estimation of the number of distinct elements of a multiset (e.g.
the distinct molecules among the reads of a sample) in a fixed amount of
memory: 4 KiB (2\ :sup:`12` one-byte registers), with a relative standard
error of about 1.6% (the estimator of Flajolet et al. 2007, with its
small-range correction by linear counting; 64-bit hash values, as in
Heule et al. 2013, so that no large-range correction is needed). The
sketch is a plain structure,
``hll_t``; sketches of parts of the data (e.g. computed by different
threads) can be merged:

.. code-block:: c

   hll_t sketch;
   hll_init(&sketch);
   for (each read)
     hll_add(&sketch, hll_hash(key, key_len));
   hll_merge(&sketch, &sketch_of_another_thread);
   double distinct = hll_estimate(&sketch);

``hll_hash()`` computes a 64-bit hash of any data; the elements are added
as hash values, so that they can be of any type.

.. _khash:

Hash tables (khash)
~~~~~~~~~~~~~~~~~~~

khash (from klib) implements hash maps and sets with macros, generating
the code for the given key and value types. It is used, for example, by
the barcode whitelist (a set of 64-bit encoded barcodes) and the cell
caller (a map from barcodes to counts):

.. code-block:: c

   #include "khash.h"

   KHASH_MAP_INIT_STR(counts, size_t)   /* string -> size_t */

   khash_t(counts) *h = kh_init(counts);
   if (h == NULL)
     fatal_printf("Failed allocating the counts\n");
   int absent;
   khiter_t k = kh_put(counts, h, name, &absent);  /* the key is not
                                                       copied */
   if (absent == -1)
     fatal_printf("Failed allocating the counts\n");
   if (absent)
     kh_value(h, k) = 0;
   kh_value(h, k)++;
   k = kh_get(counts, h, "S01");
   if (k != kh_end(h))
     printf("%zu\n", kh_value(h, k));
   kh_destroy(counts, h);

khash does not exit on allocation failures: ``kh_init()`` returns
``NULL`` and ``kh_put()`` sets its last argument to -1; they are checked
as above (or with ``LCHECK``). Other initializers:
``KHASH_MAP_INIT_INT64``, ``KHASH_SET_INIT_INT64``,
``KHASH_SET_INIT_STR``, ``KHASH_MAP_INIT_INT`` (see the header of
``khash.h``). With string keys, the strings are not copied: they must
live as long as the table (or be copied, and freed before
``kh_destroy()``).
