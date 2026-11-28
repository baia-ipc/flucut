Introduction
============

What flucut does
----------------

The whole genome of influenza is often sequenced from a multisegment
RT-PCR: universal primers bind the conserved termini of the eight
segments and amplify all of them in one reaction. For influenza A they
are two, Uni12 and Uni13; for influenza B, whose termini are shorter,
they are a cocktail of six pairs (:doc:`primers`).
The bases of a read which lie under a primer come from the primer, not
from the virus: a difference between the virus and the primer at such a
position is overwritten. flucut removes these bases from the reads,
before the assembly.

For each read, flucut:

- finds the primers at both ends of the read, in either orientation,
  with sequencing errors, also when the read begins or ends inside a
  primer;
- removes every primer-derived base and everything outside the primers
  (adapter, barcode and other flanking bases), by position;
- writes the rest of the read, in its original orientation, to a gzipped
  FASTQ file for each sample, with tags in the title telling what was
  found and removed;
- with reference sequences of the segments, tells whether each end of
  the read is at a terminus of a segment or at an interior position
  (internal priming), and reports it; it also tells, from those
  references, whether a sample is influenza A or influenza B.

Why it is needed
----------------

The untailed Uni12 primer is only twelve bases long, with a degenerate
base. Sequences similar enough to anneal to it occur inside the
segments: the primer can start an amplicon in the middle of a segment,
and the reads of such an amplicon carry the primer's bases where the
template has others. This can change the consensus at the site: a false
variant, reproducible from run to run (see :doc:`internal_priming`).
Removing the primers from the reads removes it.

IRMA, the assembler used for influenza by many laboratories, does not
remove primers: it matches a single adapter sequence, and its pipeline
MIRA does no primer trimming for influenza Nanopore reads.

What it does not do
-------------------

- It does not assign the reads to samples (demultiplexing): give it the
  reads of one sample (a file, or a directory of files, e.g. a barcode
  directory written by the basecaller).
- It does not assemble the segments or call variants: that is the work
  of IRMA or another assembler.
- It knows the primers of the multisegment RT-PCR of influenza A and of
  influenza B (:doc:`primers`); other amplicon schemes are not
  supported.
- For influenza B it trims, classifies and reports the reads as for
  influenza A, but names no subtype, reads no marker and counts no reads
  of a second subtype; its companion check against mappings,
  flucut-posthoc, is for influenza A only (:ref:`influenza-b-features`).

Platforms
---------

flucut is not specific to a sequencing platform: the errors allowed in a
primer are an option. The defaults are those of Oxford Nanopore reads.
