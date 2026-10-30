Installation
============

From the source
---------------

Only a C compiler (C11 with GNU extensions) and ``make`` are needed: the
libraries used (libcbaia, edlib, libyaml, zlib) are included in the
repository, in ``ext``.

::

   git clone https://github.com/baia-ipc/flucut
   cd flucut
   make
   make test

The program is ``bin/flucut``. ``make install`` copies it to
``/usr/local/bin`` (``make install PREFIX=<dir>`` to another prefix).

``make SYSTEM_LIBS=1`` uses the zlib and libyaml of the system instead
of the included copies.

The tests
---------

``make test`` runs the unit tests and the program tests under valgrind
(if it is installed), and the check of the coding style. The program
tests run flucut on the files of ``testdata`` and compare the results
with ``testdata/expected``.

Packages
--------

There is no Bioconda package yet.
