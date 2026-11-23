# How it works

## The idea

flucut decides from a probe: the 32 bases after a primer are located in
the references, and the read is drawn in the coverage plot from that
position and its length. That is fast, and approximate. Here each read
is aligned in full, so that every statement about a read rests on its
alignment:

- where its template begins and ends on the segment, to the base;
- whether a primer lies exactly there, and which one;
- which bases of the reference it covers, and which it lacks.

## Steps

1. **References** (`references.py`). Complete segment references are
   needed: the consensus of IRMA does not include the segment ends, and
   can be the sequence of a deleted genome or of joined amplicons. From
   a reference set of several subtypes, the reference of each segment
   is the one to which a sample of the reads maps best. It is then given
   the bases of the sample: 200,000 reads are mapped, each position
   where most reads agree on another base takes it, and this is
   repeated until nothing changes. The first and last 20 aligned bases
   of each read are left out: under a primer they are the primer's, and
   at a strong internal priming site they would write the primer into
   the reference. Insertions and deletions are not
   applied, so positions stay those of the database record.
2. **Mapping** (`mapping.py`). minimap2, with short k-mers (the
   references are a few kilobases; short templates must map), without
   secondary alignments, and without joining chains over gaps of more
   than 100 bases: the two flanks of a deleted genome are then reported
   as two alignments, and joined in the next step.
3. **Amplicons** (`amplicons.py`). The alignments of a read are grouped
   into amplicons. Consecutive alignments on the same segment and
   strand, in order, with no primer between them, are one template; if
   the reference between them has at least 100 bases more than the
   read, that is a deletion. Bases of the read which fit both flanks of
   a deletion are counted once, with the second flank. A flank of
   fewer than 30 bases cannot be mapped: where an alignment has no
   primer, a primer 12 to 60 bases before it, with the bases between
   them at the terminus of the segment, is such a flank. A primer between two alignments means two
   amplicons joined end to end.
4. **Primers** (`primers.py`). At each end of an amplicon the conserved
   terminus of each primer (Uni12 `AGCRAAAGCAGG`, Uni13
   `AGTAGAAACAAGG`) is searched where the alignment begins, with up to
   20% errors. Only the terminus is searched, not a tail, so the result
   does not depend on the primer scheme. The primer may lie inside the
   alignment (its bases match the reference at a terminus, and mostly
   at an internal site) or just before it; near a terminus, where a
   sample can differ much from the reference and the alignment begins a
   few bases inside the segment, anywhere out to the terminus. The position of the end is
   that of the first base after the primer, read from the alignment.
5. **Origin**. An end is terminal within 10 bases of the terminus of
   the segment, internal otherwise. If the reference lacks its
   terminus, bases of the read beyond the reference end are allowed for.
6. **Coverage** (`pipeline.py`). The aligned intervals of each
   amplicon, by category.
7. **Comparison** (`flucut_tags.py`, `flucut_summary.py`, `irma.py`).
   flucut's tags are compared with the amplicons of the same read. The
   category flucut gave a read is its `ct` tag (with the strand, `fs`);
   for reads written without that tag it is worked out from the tags of
   the ends. Its
   positions and its coverage (the member `coverage` of its JSON
   summary), given on its own references, are
   transferred to the references used here by aligning the two.
   IRMA's depth is transferred from its consensus in the same way.
8. **Agreement** (`agreement.py`, `origins.py`, `subtype.py`). Each
   comparison becomes a check with a share, what it is a share of, and
   limits which give its state (✓, !, ✕). The subtype of the mapping is
   the H of the HA reference chosen and the N of the NA reference
   chosen; it is compared with the subtype flucut assigned.

## Chance matches

A primer of 12 bases with two errors allowed is sometimes found by
chance. Each search is therefore repeated with control sequences: the
bases of the primer in another order, at least 5 edits from the primers
and their usual tails. Their hits are chance matches, and their rate is
reported. A control can resemble an adapter of a library; three are
used and the middle rate is reported.

## How it was checked

On reads simulated by `benchmark/simulate.py`, which records where each
base of each read comes from (14,000 reads of seven segments, 4% errors
on average, 40% of the HA reads internally primed):

| measure | result |
|---|---|
| ends with a primer at a terminus, found as such | 24,631 of 25,164 (97.9%) |
| ends with a primer at the internal site, found as such | 741 of 764 (97.0%) |
| terminal taken for internal, or internal for terminal | 0 |
| template length, both primers found: exact / within 2 bases | 88.2% / 99.8% |
| the internal priming site | at its position, 745 reads |
| primer dimers recognized | 116 of 128 |

The primers not found are those with more than two errors (the reads of
the simulation have more errors than R10 reads). Run it with
`tests/simulation_check.py`; a small version is part of the tests.

## Limits

- Reads of whole amplicons (Nanopore) are assumed. Illumina fragments
  map, but most have no primer, and paired reads are not handled.
- With only the reads written by flucut, the primers are gone: see the
  tutorial.
- A reference far from the sample in its untranslated regions makes the
  ends less exact. The read identity in the report tells how close the
  adapted reference is.
- Internal priming within 10 bases of a terminus counts as terminal.
- The reads which do not map are described (primers, groups, a typical
  sequence), not identified.
- flucut's coverage is read from its JSON summary, which gives the
  depth at up to 200 positions per segment; between them it is
  interpolated. Without the summary it is
  rebuilt from the read tags, which lack the junctions of deletions and
  the position of the reads flucut placed by a probe from inside them
  or by their length (their category is known, from the `ct` tag, but
  they are not in the coverage).
- A reference of flucut too different from the one chosen here to be
  aligned with it (e.g. the other NS allele) cannot be compared: its
  coverage and positions are left out.
- The subtype is the one named by the closest reference of the set: a
  subtype without a reference in the set cannot be found.

## Code

| module | role |
|---|---|
| `seqio.py` | FASTA and FASTQ reading |
| `references.py` | references: termini, choice, adaptation, position maps |
| `mapping.py` | minimap2 and its alignments |
| `primers.py` | primers, controls, search |
| `amplicons.py` | alignments of a read to amplicons, ends, categories |
| `analysis.py` | counting, in worker processes |
| `flucut_tags.py`, `flucut_summary.py` | what flucut said, read by read and in its JSON summary |
| `irma.py` | the results of IRMA |
| `offtarget.py` | groups of amplicons which are not influenza |
| `subtype.py` | the subtype of the mapping, of flucut and of IRMA |
| `summary.py` | the JSON summary |
| `origins.py` | the four read origins; depth and state of a segment |
| `agreement.py` | the checks of flucut against the mapping, their limits |
| `report.py` | the HTML report: the page and the measured sample |
| `report_agreement.py`, `report_subtype.py` | its first part: verdicts, checks by measure and by segment, subtype |
| `report_segments.py` | the strip of a segment |
| `report_details.py`, `report_compare.py` | the tables shown on request; IRMA |
| `plots.py`, `html.py`, `style.py` | plots, HTML helpers, stylesheet |
| `overview.py` | the table of several samples |
| `segment_coverage.py` | the page of one segment in several samples |
| `pipeline.py` | the run |
