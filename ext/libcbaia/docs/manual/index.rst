libcbaia manual
===============

libcbaia is a library of C modules for **command line programs**, in
particular **bioinformatics tools** working with sequencing reads: memory
allocation and logging with clear error messages, a uniform command line
interface, file operations, data structures, external sorting, parallel
gzip and BGZF compression, FASTQ and BAM input and output, barcode
correction, HTML reports and other building blocks of read processing
tools.

This manual is written for **C developers** who want to use the library
in their programs, or to contribute to it. The first chapters show how to
build the library, how to include it in a project and how to write a
first program; the chapters on the components describe each module, with
examples; the last chapters describe recipes for common tasks (e.g. a
demultiplexer, multithreaded processing), the tests, and how to add new
modules.

.. toctree::
   :maxdepth: 2
   :numbered:
   :caption: Contents

   introduction
   installation
   quickstart
   conventions
   core
   ds
   io
   seq
   barcodes
   recipes
   testing
   performance
   contributing
   reference
