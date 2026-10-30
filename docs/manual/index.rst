flucut manual
=============

flucut removes the primers of the influenza A multisegment RT-PCR
(the universal Uni12 and Uni13 primers) from sequencing reads, before
the assembly of the segments, e.g. with IRMA. It finds the primers at
both ends of each read, in both orientations, also when a read begins
inside a primer, and removes them by position. It recognises the reads
started by a primer at an internal site of a segment (internal
priming), whose primer bases otherwise reach the consensus as false
variants, such as the I348K change of the H5 hemagglutinin.

This manual explains how to install and use flucut, from a first run on
the test data of the repository to the classification of the reads
against reference sequences. The chapters can be read in order, as a
tutorial, or used as a reference: the last chapters list every message
and every option.

.. toctree::
   :maxdepth: 2
   :numbered:
   :caption: Contents

   introduction
   installation
   quickstart
   primers
   how_it_works
   internal_priming
   output
   reports
   live
   irma
   performance
   troubleshooting
   cli_reference
   glossary
   about
