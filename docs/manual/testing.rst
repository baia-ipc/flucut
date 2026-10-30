.. _testing:

Testing
-------

The library and the projects using it are tested at three levels: unit
tests of the modules, tests of the programs on small input files, and
further checks (sanitizers, static analysis, fuzzing).

Unit tests
~~~~~~~~~~

Each module with non-trivial logic exports a ``<module>_test()`` function
(e.g. ``whitelist_test(tmpdir)``; with a directory for temporary files,
if it needs files), using the macros of ``test.h``:

``test_assert(cond)``
   exit with an error message (file, line, condition) if the condition
   is false
``test_eq(found, exp, type, fmt)``
   compare two values of a type; output ``[OK]`` or the differing values
``test_eq_size``, ``test_eq_int``, ``test_eq_long``, ``test_eq_char``
   shortcuts of ``test_eq`` for ``size_t``, ``int``, ``long``, ``char``
``test_eq_str(found, exp)``
   compare two strings
``test_eq_mem(found, len, exp)``
   compare ``len`` characters, not 0-terminated (e.g. the content of a
   ``strbuf_t``), with a string

The unit tests of the modules reading or writing BAM files share the
helpers of ``test-bam.h`` (internal), which write BAM data or files with
a given header and reads.

e.g.:

.. code-block:: c

   void heap_test(void) {
     heap_t *h = heap_new(cmp_ints, 1);
     int v[] = {3, 1, 2};
     for (size_t i = 0; i < 3; i++)
       heap_push(h, &v[i]);
     test_eq_size(heap_size(h), 3);
     test_eq_int(*(int *)heap_pop(h), 1);
     heap_delete(h);
   }

The tests of the library are run by ``bin/cbaia-tests`` (source
``tests/cbaia-tests.c``), with one option for each group of tests, so
that a group can be run alone (``bin/cbaia-tests -w``: the whitelist);
without options, all tests are run. ``make test`` runs it under valgrind
(``--error-exitcode=1 --leak-check=full
--errors-for-leak-kinds=definite``): memory errors and leaks make the
tests fail. A new module adds its test function to the table of the test
groups, with an option letter which is not used by the common options
(``q v D L h V``), by ``-T`` or by another group; the option string and
the list of the test options in the help are generated from the table.

The projects have their own test program (e.g. ``bin/lp-tests`` of
Lumipore), with the same structure.

Program tests
~~~~~~~~~~~~~

The programs are tested on small input files, comparing their output
with the expected results. ``scripts/test-cmd-lib.sh`` provides the
``run_cmp`` bash function, for the test targets of a Makefile::

   source scripts/test-cmd-lib.sh
   run_cmp "readcount" "valgrind -q --error-exitcode=1" \
     "bin/readcount testdata/reads.fastq" 1 \
     out/readcount.tsv testdata/readcount.tsv.expected 0

The arguments: a label, a prefix of the command (e.g. valgrind), the
command, whether its standard output is redirected to the first result
file, then up to three pairs of result and expected files (the first
pair with a flag: compare after sorting). The function checks the exit
code and compares the files; on differences, the output is copied to
``<expected>.found``, for inspection (and, if correct, to become the new
expected file).

The example programs of this manual are tested in the same spirit by
``make -C docs/manual/examples check`` (part of ``make test``).

When a program reimplements an existing tool (e.g. a prototype in
Python), the expected results are produced by that tool: the C program
must produce byte-identical outputs on the same inputs.

Further checks
~~~~~~~~~~~~~~

Before a release, and after changes to the parsers or to the memory
handling:

* **Sanitizers**: build and run the tests with
  ``DEBUGFLAGS='$(DEBUGFLAGS_SANITIZE)'`` (address and undefined behaviour
  sanitizers).
* **Warnings**: the code compiles without warnings with ``-Wall
  -Wextra``, with gcc and clang (``make CPPFLAGS=-Werror``).
* **Static analysis**: e.g. ``clang-tidy
  --checks='clang-analyzer-*,bugprone-*'``, reviewing the findings (some
  are false positives, e.g. in khash).
* **Fuzzing**: see below.

Fuzzing
~~~~~~~

The parsers of the input files are fuzzed with libFuzzer (clang, with
the address and undefined behaviour sanitizers; on Debian/Ubuntu, the
clang runtime package, e.g. ``libclang-rt-18-dev``, is needed). The
harnesses are in ``tests/fuzz``:

``fuzz-reads.c``
   files of reads, through the reads reader (FASTQ, gzipped FASTQ or
   BAM): the records are read in small batches, written again as FASTQ
   and as unaligned BAM, and reverse complemented; the reads are counted
   and the completeness of the file is checked;
``fuzz-bam.c``
   BAM files: header and records, encoded again as unaligned BAM;
``fuzz-whitelist.c``
   barcode whitelists: barcodes searched at both ends with distances 0
   to 2, and cells identified from the exact matches.

::

   make fuzz                         # builds bin/fuzz/fuzz-<harness>
   tests/fuzz/run.sh reads 600 4     # harness, seconds, parallel jobs

The corpus is kept in ``obj/fuzz/corpus-<harness>`` (started from the
seeds in ``tests/fuzz/seeds/<harness>``); crashing inputs are written to
``obj/fuzz/crashes-<harness>/`` and can be replayed with
``bin/fuzz/fuzz-<harness> <file>``.

As the parsers exit the program on invalid input, ``exit()`` is wrapped
(linker option ``--wrap=exit``) so that it returns to the harness
(``FUZZ_RUN()`` of ``tests/fuzz/fuzz-common.h``), which closes the files
left open; leak detection is disabled, as the memory of an interrupted
parse is not freed. ``fuzz-common.h`` is also used by the harnesses of
the projects, which fuzz their own parsers (e.g. the sample sheets and
library configurations of Lumipore): they add ``tests/fuzz`` of the
library to the include path, and define ``FUZZ_NAME`` (used in the name
of the temporary input file).
