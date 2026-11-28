# Reading the report

`posthoc_report.html` is one self-contained file (inline plots, light and
dark themes). It goes from the verdicts to the details: read the top,
and open the rest when you need it. The page explains its own terms in
"Terms used in this page", at its end.

Each section and each box folds: click its title (▾ open, ▸ folded). The
sections open at the start, except "All numbers"; a link to a part of
the page opens the folded parts around it.

Conventions:

- **measured** is what the mapping shows; **flucut** (or predicted) is
  what flucut reports without mapping.
- Read origins have the colours of the report of flucut: full template
  (blue), internally primed (magenta), internal deletion (violet),
  other (grey).
- A state always has a symbol: `✓` within the limit, `!` within the
  wider limit, `✕` outside, `·` not judged (fewer than 20 cases, or
  nothing to compare).
- Every percentage names what it is a share of. `0.0%` and `100.0%`
  mean exactly none and exactly all: a share which is not zero, but
  below 0.1%, is written `<0.1%`, and one above 99.9%, but not all,
  `>99.9%`.

## flucut and the mapping

Shown when the reads or the summary of flucut are given.

**The verdicts.** Two lines of badges:

- *Subtype*: `✓ Subtype H5N1: confirmed by the mapping`,
  `✕ Subtype: flucut H5N1, mapping H5N8` or `· Subtype: cannot be told`,
  followed by the subtype of the mapping and of each IRMA run.
- *Checks*: how many checks are outside their limits, and which.

**Tiles**: the main checks at a glance.

**Subtype.** The subtype of the mapping is the H of the HA reference
the reads map to best and the N of the NA reference. For HA and NA the
table gives the type and the reference of the mapping, the reads behind
it (of the first 20,000 reads, those which map to references of that
type), the next type with its reads (a mixed sample, or a weak call),
the type and the reference of flucut, and the type IRMA assembled. With
fewer than 20 reads on HA or NA the subtype cannot be told. The other
segments are listed with the two references and the subtype of their
strains, as information: another reference for an internal gene is not
another subtype.

**Agreement by measure.** One line per check, with its limits:

| check | meaning |
|---|---|
| same origin | of the ends both call terminal or internal, those called the same |
| internal / terminal calls confirmed | of the ends flucut calls so, those the mapping places there |
| internally primed / terminal ends found | of the ends the mapping places there, those flucut calls so |
| position within 20 bases | of the ends located by both |
| primer bases left | of the internally primed ends flucut trimmed: these bases reach the consensus |
| same read origin | of the amplicons compared read by read, those in the same of the four origins |
| same number of amplicons in the read | of the reads written by flucut which map |
| depth, flucut / measured; difference per position | all segments together |

**Agreement by segment.** State (complete, partial, low depth, not
found), median and mean depth, measured and of flucut; then the ratio
of the two depths, their mean difference per position and their
correlation along the segment. A segment agrees when the state is the
same and the depth is within its limits.

## The sample, measured

- **Badges and tiles**: segments complete, what deserves attention
  (amplicons of another template, internal priming, internal
  deletions), and the main numbers.
- **Reads**: what the reads are (influenza, other amplicons, not
  mapped, primer dimers), the amplicons by origin and the primers
  found, each next to the same bar of flucut.
- **Internal priming sites**: the largest sites, with the reads flucut
  reports at the same place, and the template under the primer. A site
  whose template differs from the primer at more than 5 bases is not
  one where the primer annealed: the primer is joined to the segment
  there, often at the junction of a deleted genome whose other flank is
  too short to align.
- **Reads assembled by IRMA** (with `-i`): for each kind of read, the
  share each IRMA run assembled.
- **Segments**: one strip per segment. At the left, the numbers
  measured next to those of flucut, and one line for each consensus of
  IRMA. In the middle, the coverage by read origin measured and, below
  it on the same scale, the one flucut predicts; the lines are the
  depths of IRMA, the vertical line the largest internal priming site,
  the red bar the positions below the depth of a complete segment. At
  the right, the lengths of the complete amplicons (bars: measured;
  line: flucut).

## All numbers

Folded at the start. Every table, folded too: references, ends of the amplicons, primer search
(errors, chance matches measured with control sequences, bases before
the primers), internal priming sites with the template under the
primer, complete amplicons, amplicons by category, internal deletions,
ends without a primer, amplicons which are not influenza; then the
comparison with flucut read by read (origin and position of the ends,
bases removed, counts of the run, amplicons in a read, read origin) and
the reads of each IRMA run by step and by kind.

## Overview of several samples

`overview.html` (`flucut-posthoc-overview`) has one line per sample:
the verdict, the subtype (badge, flucut, mapping, IRMA), the segments
which agree, the measures and the checks, with cells tinted by their
limits. `overview.tsv` has the same values; its last columns are
`subtype_flucut`, `subtype_mapping`, `subtype_irma` and
`subtype_status` (`same`, `different`, `unknown`).

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
by category, flucut's, and IRMA's), `irma`, `flucut`, `offtarget` and
`subtype`: its `status` (`same`, `different`, `unknown`) and `text` (the
verdict), the subtypes `measured`, `flucut` and `irma`, for HA and NA
(`genes`) the type, the reads by type and the next type, and for each
segment (`segments`) the reference of the mapping and of flucut.
