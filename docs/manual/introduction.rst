.. _introduction:

Introduction
------------

Why a library
~~~~~~~~~~~~~

Command line tools written in C, in particular in bioinformatics, need
again and again the same basic functionality: allocating memory and
exiting with a clear message if this fails, logging messages with
different levels of verbosity, opening and reading files with
descriptive error messages, parsing the command line options in a
uniform way, writing unit tests, reading FASTQ and BAM files, writing
compressed output in parallel, and some recurring data structures and
algorithms.

Instead of copying these modules from one project to the next (and
fixing the same bugs in each copy), libcbaia collects them in a single
place. It is used by the C projects of the bioinformatics group (BAIA)
of the Institut Pasteur du Cambodge, e.g. Lumipore, a demultiplexer of
Illumina libraries sequenced on Nanopore: most of the generic code of
Lumipore (FASTQ and BAM input, output files, barcode correction, cell
calling, library complexity) is part of libcbaia, so that other tools
can use it.

Design
~~~~~~

The library follows a few principles, which make the programs using it
short and robust:

**Errors exit the program, with a clear message.** Memory allocation,
file operations and the parsers of the input formats do not return
error codes: if something goes wrong (a file cannot be opened, the
memory is exhausted, the input is invalid), the program is interrupted
with a message describing the problem in terms of the user's input
(file name, line, value). The caller thus does not need to check the
return values, and the code stays focused on the task.

**Messages go through a logger.** The modules never write to the
standard error directly: they receive a ``Logger``, created by the
program from the command line options (quiet, verbose, debug, log file),
and use it for errors, warnings and information. The standard output is
reserved for the results of the program.

**Small, independent modules.** Each module encapsulates one data
structure or one task, behind a small interface: an opaque type with a
constructor (``_new()``), functions operating on it and a destructor
(``_delete()``). The modules have no global state, so that they can be
used by several threads, each with its own objects.

**Few dependencies.** The library only needs the C standard library,
POSIX and zlib (included in the repository, so that the library can be
built with only a C compiler and ``make``). In particular, BAM files are
read and written without htslib.

**Tested.** Each module has unit tests, run under valgrind (failing on
memory errors and leaks); the parsers of the input files (FASTQ, BAM,
whitelists) are fuzzed.

Components
~~~~~~~~~~

The modules are organized in four components (directories of ``src``),
shown in :numref:`fig-components`:

``core``
   the base functionality of a command line program, used by all other
   components: memory allocation (``alloc.h``), logging (``logger.h``),
   command line interface (``argparse-def.h``), files and directories
   (``file.h``), strings (``strutils.h``, ``strbuf.h``), numeric
   utilities (``numutils.h``), progress display (``progress.h``),
   terminal colors (``colors.h``), unit tests (``test.h``);

``ds``
   data structures: bit operations, Bloom filter, set of integers, heap,
   HyperLogLog;

``io``
   handling of large files: sorting of files which do not fit in memory,
   parallel gzip and BGZF compression, MD5 checksums; self-contained HTML
   reports with plots;

``seq``
   sequences and sequencing reads: reverse complement, FASTQ and BAM
   reading and writing, a common reader of FASTQ and BAM files, SAM tags
   and MinKNOW metadata in the read titles, the header of BAM output,
   output files of several destinations, the edit distance of short
   sequences, the k-mer index of a set of sequences, barcode whitelists,
   cell calling and library complexity.

.. _fig-components:

.. figure:: figures/components.*
   :width: 100%

   The components and modules of libcbaia. The arrows show the main
   dependencies between modules (all modules use ``alloc.h`` and
   ``logger.h`` of ``core``); external code is shown in grey.

The header names are unique over all components: the modules include
each other, and are included by the programs, by their name (e.g.
``#include "logger.h"``), with the directories of the components in the
include path.

What is not in the library
~~~~~~~~~~~~~~~~~~~~~~~~~~

The library contains generic functionality, which is useful to several
programs. Code specific to one tool (e.g. the analysis of the structure
of Illumina library molecules in Lumipore) stays in the tool. Code which
needs further libraries (e.g. alignment with edlib, YAML parsing with
libyaml) also stays in the tools which use it, so that the library keeps
its minimal dependencies. :ref:`contributing` describes when and how to
move code into the library.
