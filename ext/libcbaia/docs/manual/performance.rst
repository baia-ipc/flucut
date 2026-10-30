.. _performance:

Performance
-----------

Read processing programs handle millions of records; the library is
designed so that the programs using it spend their time in their own
work, not in memory management, input or output.

Guidelines
~~~~~~~~~~

**No allocations in the inner loops.** The batches of records, the
output buffers (``strbuf_t``) and the compressors are allocated once
(per thread) and reused: they keep their memory when cleared, and grow
only when a larger record arrives. The memory usage does not depend on
the size of the input (e.g. Lumipore uses 20–50 MB, whatever the size
of the run).

**Batches.** The readers return batches of records stored in one
buffer, which amortizes the cost of the input and of the locking in a
multithreaded program (a batch of 1,000–10,000 records is a good size).

**Parallel compression.** Writing gzipped output is often the most
expensive step: with the gzip and BGZF compressors, each thread
compresses its own output (:ref:`recipe-threads`).

**Cache-friendly structures.** Arrays with ``len`` and ``cap`` fields
are preferred to linked structures; barcodes are encoded in 64-bit
integers (2 bits per base) and looked up in hash sets; the HyperLogLog
sketches count distinct elements in 4 KiB.

**No logging in the inner loops.** Messages are logged per file or per
phase, not per record; ``TRACE`` and ``DEBUG_TRACE`` are for temporary
debugging only.

**Measure before optimizing.** Use realistic inputs, and compare the
versions with the same data and machine (alternating the runs, and
taking the median, as the time varies between runs), or count the
instructions with ``valgrind --tool=callgrind``, which does not depend
on the load of the machine. For profiling, build with
``DEBUGFLAGS='$(DEBUGFLAGS_PROF)'`` (gprof) or use ``perf``.

Link-time optimization
~~~~~~~~~~~~~~~~~~~~~~

As the library is a static library, its functions are not inlined into
the program by default. Compiling the program and the library with
``-flto`` (as the ``DEBUGFLAGS_NDEBUG`` of Lumipore does) allows the
compiler to inline small functions of the library (e.g. those called
for each read) into the loops of the program, as if they were in the
same file. The vendored library is built with the ``DEBUGFLAGS`` of the
project, and rebuilt automatically when they change (e.g. after a build
of the library alone with other flags, by ``make test`` in its
directory).
