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
                             (``-e auto``), the worker threads and the
                             live mode
``src/seq``                  primers, the scheme, references and
                             segments; the search of primers and
                             adapters; the analysis of an amplicon
                             (chimeras, internal priming, junctions,
                             placement) and the trimming of the reads
``src/analysis``             the counters and what is computed from
                             them: lengths, coverage by read origin,
                             quality summary and flags, subtype reads,
                             the cleavage site
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

As in libcbaia, a header is included by its name, without the
directory (``#include "counter.h"``): the names are unique over all
the directories, which the Makefile passes to the compiler, and the
objects are all in ``obj``. To add a module, put it in the directory of
its purpose; the Makefile finds it (the list of directories is
``COMPONENTS``). ``tests/check-style.sh``, run by ``make test``,
checks all of them.

Packages
--------

There is no Bioconda package yet.
