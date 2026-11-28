Installation
============

From the source
---------------

Only a C compiler (C11 with GNU extensions) and ``make`` are needed: the
libraries used (libcbaia, edlib, libyaml, zlib, libdeflate) are
included in the repository, in ``ext`` (their versions and sources:
``ext/README.md``).

::

   git clone https://github.com/baia-ipc/flucut
   cd flucut
   make
   make test

The program is ``bin/flucut``. ``make install`` copies it to
``/usr/local/bin`` (``make install PREFIX=<dir>`` to another prefix).

``make SYSTEM_LIBS=1`` uses the zlib and libyaml of the system instead
of the included copies. libdeflate, which compresses the output, is
always the included copy.

The tests
---------

``make test`` runs the unit tests and the program tests under valgrind
(if it is installed), and the check of the coding style. The program
tests run flucut on the files of ``testdata`` and compare the results
with ``testdata/expected``.

The source code
---------------

The C sources are in ``src``. The programs are at its top level
(``flucut.c``, the command; ``fc-tests.c``, the unit tests); the modules
are in subdirectories, by purpose:

===========================  ===========================================
Directory                    Content
===========================  ===========================================
``src/core``                 the options and their parsing, the banner,
                             the settings decided at the start of a run
                             (``-e auto``), the worker threads, the
                             trimming of a batch of reads and the live
                             mode
``src/seq/primers``          what is searched in a read: the built-in
                             primers, the primer schemes, the names of
                             the segments and the types of influenza
``src/seq/references``       the reference segments: the built-in set,
                             the subtype of a reference (from its FASTA
                             description), the k-mer index and the
                             location of a probe
``src/seq/search``           how a pattern is searched in a read: the
                             alignment of a primer (whole, at an end of
                             the read, at a segment terminus), of the
                             Nanopore adapters and of an anchor site
``src/seq/amplicon``         the analysis of the primers of a read: the
                             primers at its ends, the primer matches
                             inside it, the origin of each end and the
                             part to keep
``src/seq/reads``            the rest of the structure of a read: its
                             placement on the references, the junction
                             of an internal deletion and the split of a
                             chimera
``src/seq/markers``          the sites read from inside a read: the
                             table of the markers, their anchors in the
                             references and the reading of one read
``src/analysis``             the counters and what is computed from
                             them: lengths, coverage by read origin,
                             quality summary and flags, subtype reads,
                             the cleavage site and the markers
``src/io``                   the groups of input files, the output
                             files of the reads, the state saved for
                             ``--resume`` and the result tables
``src/reports``              the content of the reports and their
                             writing, shared by the formats
``src/reports/html``         the HTML report: its parts, the templates
                             (``*.html``, embedded in the program) and
                             the logo data
``src/reports/json``         the JSON summary
===========================  ===========================================

In ``src/seq``, the references and the searches use the primers;
``reads`` and ``markers`` use those three, and ``amplicon`` all four.
No directory uses one which depends on it.

The line between ``src/seq`` and ``src/analysis`` is the read: a module
of ``src/seq`` works on one read, or on one reference, and keeps
nothing between them; a module of ``src/analysis`` counts what the
modules of ``src/seq`` found, over the reads of a sample, and says what
the counts mean. A site therefore has its anchors and its reading in
``src/seq/markers`` and its counts and its verdict in
``src/analysis``.

As in libcbaia, a header is included by its name, without the
directory (``#include "counter.h"``): the names are unique over all
the directories, which the Makefile passes to the compiler, and the
objects are all in ``obj``. To add a module, put it in the directory of
its purpose; the Makefile finds it (the list of directories is
``COMPONENTS``). ``tests/check-style.sh``, run by ``make test``,
checks all of them.

A large module is split into files by step, with one header for the
rest of the program. The analysis of a read is used through
``amplicon.h``; its steps (``amplicon-ends.c``, ``amplicon-cut.c``,
``amplicon-adapters.c``, ``amplicon-locate.c``) share the state of the
analyzer through ``amplicon-analyzer.h``, which only they include. The
unit tests of such a module are in a file of their own
(``amplicon-test.c``, ``reference-test.c``).

Packages
--------

There is no Bioconda package yet.
