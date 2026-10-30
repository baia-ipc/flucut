# flucut

Primer trimming for influenza A whole-genome sequencing.

flucut removes the primer-derived bases from the reads of the influenza A
multisegment RT-PCR with the universal Uni12/Uni13 primers, before the
assembly of the segments (e.g. with IRMA). It is written for Oxford
Nanopore reads; support for Illumina reads is planned.

## Why

The universal primers bind the conserved termini of the segments, which
are the same in all eight segments and in all subtypes. The bases of a
read which lie under a primer therefore come from the primer and not
from the template: a difference between the template and the primer at
such a position is overwritten, and a variant seen there is an artefact.

The Uni12 primer is short: used without a tail it is only the twelve
bases of the conserved terminus, with a degenerate position. A sequence
that similar also occurs inside the segments, where the primer can
anneal and start an amplicon in the middle of a segment. The reads of
such an amplicon begin with primer-derived bases which are not the
template there, and the primer writes its own sequence into the
consensus: a variant which looks real and is not, reproducible from run
to run, and present in the public databases. This is what flucut is
built around; see `docs/internal-priming.md`.

IRMA does not trim primers (it matches a single adapter sequence,
exactly or with one mismatch), so the trimming has to be done before it.

## What it does

- finds the Uni12/Uni13 primers at both ends of a read and in both
  orientations, whole or partial (a read can begin inside a primer),
  with the error rates of Nanopore reads;
- removes, by position and not by difference, the primer-derived bases;
- because flucut knows which primers were used, it knows where inside a
  primer the tail added by the protocol ends and the conserved terminus
  begins, and can therefore offer a choice no generic trimmer can:

  - remove every primer-derived base, tail and conserved terminus
    together (the default): no base of the output comes from a primer,
    and the consensus is shorter than the segment by the termini
    (12 and 13 bases);
  - remove only the tail and keep the conserved terminus: the consensus
    has complete termini, but those bases come from the primer.

- a primer is accepted as the primer of an end only near that end of
  the read: a match deeper inside is not the primer which started the
  amplicon, and trimming there would remove template. Those matches are
  reported instead, as the sites at which the primer can anneal and as
  the junctions of chimeras;
- reads which are not amplicons are counted and reported: primer dimers,
  reads without a primer, reads with a primer at one end only.

## The primers

The conserved termini and the primers built into flucut are those of the
publications, quoted and not derived; each sequence is documented with
its source in `src/primers.h`. The built-in schemes are the tailed
primers of the multisegment RT-PCR (`mbt`, the default) and the
conserved termini used as primers without a tail (`uni`), and the
combination used in practice of an untailed degenerate Uni12 primer with
a tailed Uni13 primer (`uni12-mbt13`, the default).

## Building

Only a C compiler and `make` are needed: the libraries used (libcbaia,
edlib, libyaml, zlib) are included in `ext`.

```
make
make test
```

`make SYSTEM_LIBS=1` uses the zlib and libyaml of the system instead of
the included ones.

## State

Early development. The parts which are implemented and tested are the
built-in primers, the search of a primer in a read (with degenerate
bases and partial primers at the read ends), the primer schemes, the
analysis and trimming of a read, and the detection of the primer
matches inside a read. The command line program, the reports
and the live mode are being written.

## License

ISC (see `LICENSE`).
