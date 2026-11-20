# Reading the report

`posthoc_report.html` is one self-contained file (inline plots, light and
dark themes). In it, `0.0%` and `100.0%` mean exactly none and exactly
all: a share which is not zero, but below 0.1%, is written `<0.1%`, and one
above 99.9%, but not all, `>99.9%` (`<0.01%` and `>99.99%` where two
decimals are shown). Its parts, in order:

## What the reads are

Every read is mapped. A read which does not map is searched for the
primers at its two ends:

| class | meaning |
|---|---|
| influenza amplicons | at least one part of the read maps to a segment |
| other amplicons | both primers, a template between them, no mapping: an amplicon of another template |
| not mapped, one primer / no primer | the rest of the reads which do not map |
| primer dimers | both primers and fewer than 30 bases between them |

A read can contain several amplicons joined end to end; each is analysed
on its own ("amplicon" below), as flucut does when it splits a read.

## References

The reference used for each segment, how many of its bases were replaced
by those of the sample, the identity of the reads to it, and, for each
IRMA run, how its consensus compares: its length, the share of the
reference it covers, its deletions, and primers found inside it.

## Ends of the amplicons

Left and right are those of the segment (the start and the end of its
plus strand), whatever the orientation in which a read was sequenced.

- **terminal**: the template begins (or ends) within 10 bases of the
  terminus of the segment (`--terminal-window`);
- **internal**: anywhere else;
- with a **primer** before it, an internal end is internal priming;
  without one, it is the end of a fragment.

"How well the primers were found" gives the errors in the primers found
and the rate of chance matches, measured with control sequences. "What
precedes the primers" shows the bases before each primer: this is where
you see which tailed primer, barcode or adapter the library really has.

## Internal priming sites

The internal ends with a primer, grouped into sites as in flucut
(positions at most 5 bases apart). For each site: the sequence of the
reference under the primer and its mismatches with the primer, and the
reads flucut reports at the same place.

## Complete amplicons

The amplicons with Uni12 at the left and Uni13 at the right, by segment:
their number, the usual length of the template, those shorter than 90% of
the reference. The histograms have the lengths measured (bars) and those
of flucut (line).

## Coverage by read origin

The depth along each segment, stacked by what the amplicons are:

| band | amplicon |
|---|---|
| terminal + / − | a primer at the left end, at the terminus; read on the plus / minus strand |
| internally primed + / − | a primer at the left end, inside the segment |
| internally primed at the 3′ end | the left end terminal, a primer at the right end inside the segment |
| no primer at the 5′ end | starts at the terminus, or inside the segment, without a primer |
| internal deletion | at least 100 bases of the reference missing (`--min-deletion`); only the bases present are covered |

The lines are the depths reported by IRMA. At the right of each plot is
the plot of flucut for the same segment, and below them a table: the
mean depth predicted and measured, their correlation along the segment,
and the mean difference per position.

## flucut compared with the mapping

- **Origin of the ends**: flucut's call for each end (terminal,
  internal, unknown, no primer) against the mapping. "Confirmed" is the
  share of flucut's calls which the mapping agrees with; "found by
  flucut" the share of the ends of the mapping which flucut calls so.
- **Position of the ends**: how far flucut's position is from that of
  the mapping. A reference of flucut which lacks bases at its start
  cannot give positions before its first base.
- **Bases removed**: at each end with a primer, the bases flucut removed
  against the bases before the template. *Primer bases left* at an
  internally primed end is the error which matters: those bases reach
  the consensus.
- **Counts of the run**, **amplicons in a read**, **category in the
  coverage plot**: the same comparison for the counts of the summary,
  the splitting of reads, and the band a read is drawn in.

## IRMA

For each IRMA run: the reads at each of its steps, and, for each kind of
read, how many IRMA assembled.

## Other tables

The most frequent internal deletions, the positions where many amplicons
end without a primer, and the amplicons which are not influenza.

## Output files

`reads.tsv.gz` has one line per amplicon of each read (and one line for a
read which does not map):

| column | content |
|---|---|
| `read`, `unit`, `units` | the read, the number of the amplicon in it, the amplicons in it (0: not mapped) |
| `read_length` | bases of the read |
| `reference`, `strand` | the segment reference and the strand of the read |
| `left_position`, `right_position` | first and last template base on the reference (1-based) |
| `left_primer`, `right_primer` | the primer before the left end, after the right end (`.`: none) |
| `left_primer_errors`, `right_primer_errors` | errors in the primer |
| `left_origin`, `right_origin` | `terminal` or `internal` |
| `template_length` | bases of the read between the two ends |
| `deletions` | deleted intervals of the reference (1-based, inclusive) |
| `category` | the band in the coverage plot; for a read not mapped, its class |
| `identity` | of the aligned bases to the reference |

`posthoc_summary.json` contains everything the report shows: `reads`,
`amplicons`, `ends`, `primer_search`, `internal_sites`,
`complete_amplicons`, `coverage` (for each reference, the measured depth
by category, flucut's, and IRMA's), `irma`, `flucut` and `offtarget`.
