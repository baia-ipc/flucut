# Internal priming

This is the problem flucut is built around. It is described here in one
place, because it drives the default settings, the reports and the
benchmark.

## The artefact

The universal primers of the influenza A multisegment RT-PCR bind the
conserved termini of the segments. The Uni12 primer, in the form used
without a tail, is only the twelve bases of the conserved terminus, with
a degenerate position:

    AGCRAAAGCAGG

Twelve bases, one of them degenerate, are short enough to occur, or
nearly occur, by chance inside a segment of one or two thousand bases.
The primer then anneals there and starts an amplicon in the middle of
the segment. The reads of that amplicon begin with the bases of the
primer; at that position those bases are not the template. They are read
as if they were sequence, and if enough reads carry them, they become
the consensus.

The result is a variant which looks real and is not. It is reproducible,
because the site is the same in every run, and it survives into the
public databases.

### The case which motivated flucut

In the H5N1 sequences of the Virology Unit of the Institut Pasteur du
Cambodge, the HA segment carried an isoleucine to lysine change at amino
acid 348: the motif at positions 348-354 read

    GLFGAKA   instead of   GLFGAIA

at nucleotides 1057-1060, `AAA` instead of `ATA`. Of 787 HA sequences
checked, 140 carried the wrong form.

The cause is a region inside HA which is similar enough to the Uni12
primer for it to anneal there. The primer carries its own `AAA` into the
reads, overwriting the `ATA` of the template; a single base changes, and
the codon changes from isoleucine to lysine.

The Uni13 primer, which is used with its tail and is 23 bases long, does
not do this: it is long enough to be specific.

## What flucut does about it

### Trimming removes the artefact

A read started at an internal site looks exactly like a read of a real
amplicon: the primer-derived bases are at the end of the read. flucut
trims by position at the ends of the read, so those bases are removed
from an internally primed read just as from a real one. The wrong bases
never reach the assembler.

This is why the default is `--trim primer`, which removes every
primer-derived base. **`--trim terminus`, which keeps the bases of the
conserved terminus, keeps exactly the wrong bases at an internal site
and reintroduces the artefact.** It exists for protocols where internal
priming is excluded and complete termini are needed; flucut warns when
it is used.

### A primer is only searched near the ends of a read

A match deep inside a read is not the primer which started the amplicon.
Trimming there would remove real template, which is the opposite error.
The primer of an end is therefore accepted only within a window of that
end of the read (`AMPLICON_DEFAULT_WINDOW`, overridable).

The unit tests show both sides of this: with the window, a read whose
only primer match is deep inside it is not trimmed at all; without it,
the same read loses everything before that match.

### Detection and reporting

Removing the artefact is not enough: the point is also to show that it
happened, and how often. flucut reports it in three layers; the first is
implemented, the other two are planned with the counting pipeline.

1. **Internal matches** (implemented). The matches of a primer between
   the two end primers are collected (`amplicon_t.internal`), with their
   position, which primer matched and how many errors. These are the
   sites at which the primer can anneal, and the junctions of chimeras
   and of reads containing several amplicons. They are counted, never
   trimmed. A run in which the Uni12 primer matches inside many reads is
   a run in which internal priming is likely.

2. **Start contexts** (planned). After trimming, the bases which follow
   the removed 5' primer are the first bases of whatever the amplicon
   started at. Reads of real amplicons start at one of sixteen places
   (eight segments, two ends). Internally primed reads start somewhere
   else, and always at the same place, so they form their own group.
   Grouping the reads of a run by the first k bases after the trimmed
   primer, and ranking the groups, therefore separates the real starts
   from the internal ones without needing any reference sequence. The
   report names the unexpected groups and gives their sequence and their
   share of the reads, which is the number to put in a paper.

3. **Against references** (planned, optional). With the segment
   references of the run, each read start can be classified as terminal
   or internal outright. This is the most direct measure, and the one to
   use when references are at hand.

## Consequences for the benchmark

The measure which matters is not only how many primer bases are left in
the reads, but whether the artefact reaches the consensus. The benchmark
therefore reports, for each tool:

- the primer bases left at the ends of the reads;
- the template bases wrongly removed;
- **the bases of the consensus which come from a primer**, at the
  termini and at internal priming sites, with the IRMA consensus of the
  trimmed reads compared with the truth;
- the reads lost.

A tool which trims only what it recognises as an adapter at the very
start of a read, or which trims by sequence difference instead of by
position, can leave the internal artefact in place while looking
perfect on the terminal measures. That is the comparison which decides
whether flucut is worth using.
