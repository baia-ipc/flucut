# The markers, read from the reads

flucut reports what the reads carry at a few positions of the influenza
proteins whose residue is worth knowing before anything is assembled:
the residues that decide whether a virus binds the receptors of birds or
of humans, replicates in a mammalian host, or resists a drug.

It reports them as a distribution, not as a call. For each marker it
gives the residue, the reads it was read from, and the share of those
reads carrying each value. A residue present in a part of the reads is
therefore visible as such — which is the difference from reading the
same position off a consensus, where a minority is rounded away.

Nothing here is assembled or mapped. Each position is found by a short
conserved peptide directly beside it, searched in the read itself
(`src/seq/search/anchor-site.h`), so the markers are available as the
reads arrive.

## Adding or changing a marker

The table is not in the program sources. It is `metadata/markers.yaml`, one
record per marker, which the build turns into a header
(`scripts/markers-table.awk`). Adding a marker is a record in that file:
its anchor peptide, where the window lies beside it, the positions and
the convention they are counted in, and one line saying what the residue
means. Nothing in the C sources changes, and the column of
`summary.tsv`, the JSON object and the row of the report follow from the
record.

The one thing to check by hand is the anchor: run with `--debug` and
read the lines which name the references that do not give the new
marker.

## What is in the table, and why

A marker is in the table only if it fulfils both of these:

- **its evidence is settled**: the residue has a published phenotype
  which does not change with the naming of the clades, so that the table
  needs no upkeep;
- **it can be found locally**: a short conserved peptide lies directly
  beside it, so that the position is reached from the sequence itself.

Clades and genotypes are therefore not here, and will not be: they are
defined by nomenclatures which are revised, and a table compiled into a
program cannot follow them. A program which named a clade would be
wrong a year later without anything in it having changed.

| marker | protein | positions | numbering | what it means |
|---|---|---|---|---|
| `PB2-627` | PB2 | 627 | PB2, from the initiator methionine | lysine replicates in mammalian cells at their temperature; glutamate is the avian residue |
| `PB2-701` | PB2 | 701 | PB2, from the initiator methionine | asparagine adapts the polymerase to mammals, as 627K does and without it |
| `PA-38` | PA | 38 | PA, from the initiator methionine | threonine lowers the susceptibility to baloxavir |
| `NA-275` | NA | 275 | N1 numbering (274 in N2 numbering) | tyrosine lowers the susceptibility to oseltamivir in the N1 neuraminidases |
| `HA-226-228` | HA | 226, 227, 228 | H3 numbering of the mature HA1 | leucine at 226 with serine at 228 binds the receptors of the human airway; glutamine with glycine binds the avian ones |
| `M2-31` | M2 | 31 | M2, from the initiator methionine | asparagine resists amantadine and rimantadine, and is the residue of most viruses now circulating |
| `NA-stalk` | NA | the stalk, between residues 36 and 90 | bases from the start of NA to its conserved HxNxT peptide | a shorter stretch is a deletion in the stalk, which poultry-adapted viruses carry |
| `NS1-tail` | NS1 | the C-terminus | the residues before the stop codon | the last four residues are the PDZ-binding motif, and a stop before them is the truncated C-terminus of the viruses circulating in humans |
| `PB1-F2` | PB1-F2 | the whole reading frame | codons of PB1-F2 | a frame of 87 to 90 codons encodes the whole protein; a stop codon early in it leaves a fragment |

The markers are of influenza A. None of the anchors is found in the
influenza B references, so an influenza B sample reports no marker
rather than a wrong one.

## Numbering

The position of a residue depends on what it is counted on. The
haemagglutinin is numbered on its precursor HA0, on the mature protein,
or on the structure of H3; the neuraminidase of N1 and of N2 differ by
one at the oseltamivir position; a virus whose stalk is shorter has
every residue after it shifted.

No marker here uses a position in the sequence. The residue reported is
the one beside the anchor, whatever its number in that sample, and the
number in the table says only which published residue it is, together
with its convention. This is not a detail: in the reference set flucut
ships, the residue equivalent to H275 of N1 sits at position 272, 274 or
275 of the neuraminidase of the reference, depending on the length of
its stalk, and the residue equivalent to 226 of H3 sits at 242 in one
H3 reference and at 243 in another, because the signal peptide of the
second is one residue longer. An absolute position would read the wrong
residue in both.

## How a marker is found

The anchor of a marker is the bases of the reference which encode a
short conserved peptide directly beside the position. It is searched in
each read with the aligner used for the primers (edlib, infix, the IUPAC
codes equal to their bases), in both orientations, with the same errors
allowed as for a probe, and the window beside the match is read as the
read carries it.

The anchor is directly beside the window, so an insertion or a deletion
between the two cannot shift what is read out of frame. A window before
the anchor is counted from the first base of the match, one after it
from the last, so an indel inside the anchor moves neither.

The peptide is a pattern: a residue, a `.` which any residue matches, or
a class such as `[AS]`. The classes keep an anchor which has to cover
several subtypes specific enough not to be found elsewhere in the
segment. The anchor of the receptor-binding site, for instance, is

    R[IM][DENST][YFIHV][YFSH]W[TSALGMI][LIVM]

which is the conserved arginine of the 220-loop, four residues which
vary with the subtype, and the tryptophan after them: the three residues
before it are 226, 227 and 228 in H3 numbering. Written loosely enough
to cover the subtypes without the classes, the same pattern also matches
sequence which is not the site.

An anchor is built for a reference only if

- the peptide occurs exactly once in the reference, in one reading
  frame, and
- the bases of the anchor are found nowhere else in the reference within
  the errors a read is searched with.

Both are checked when a run starts, for every reference, and a reference
which fails either is reported at `--debug` and simply does not give
that marker. An anchor which another stretch of the segment could answer
for would report the residue of the wrong place.

## Which reads are read

A read of a whole amplicon (Nanopore) is located on a reference, so only
the markers of that reference are looked for in it, and each anchor
first near where the read's place on the reference puts it (the whole
read is searched only if it is not found there). A Nanopore read placed
on no segment gives no marker, as it gives no cleavage site: what it is
a read of is not known. A read of a fragment (Illumina) which covers a
marker usually carries no primer and is not located at all: every marker
is looked for, and the anchors themselves say which ones the read
covers.

**A marker is read only from the part of a read flucut keeps**, the
bases it judges to be the template. The bases it removes are the whole
primer at an internally primed end, the flanks and the adapters, and
reading a window over them reads the primer: the conserved terminus
`AGCAAAAGCAGG` ends in `AAA`, which is a lysine, so every marker whose
window is the codon before its anchor reads K on a read primed just
before that anchor. In one sample this turned PB2 627 into K in 89% of
423 reads where the kept bases of the same sample say E in 78% of 54.

**A read with an internal deletion gives a marker only where the anchor
and the window lie on the same side of its junction.** Such a read is
two flanks of the segment joined, and the bases before the junction come
from somewhere else entirely. A read whose junction was not located
gives no marker at all. The same
holds for the short deletions beside the junction, which are missing
stretches of the segment just as the junction is, and for a primer
which matched inside a read.

Together these say one thing: **a marker is read only from a stretch
of a read which follows its reference**. It is not enough for the
anchor to be found; the bases beside it have to be the bases of that
position. The anchor must also lie in the stretch of the reference the
read covers, so that an anchor found where the read does not reach
gives nothing.

What a marker says does not depend on how the reads are written: the
cuts and the kept part are those of the read as it arrived, whether or
not `--deletions` writes it in pieces, and a read a category filter
leaves out of the output is still read for its markers. Which bases of
a read are the template is a property of the data; which reads are
written is a choice.

## How far it can be trusted

Three numbers are reported for every marker, and a verdict made from
them.

- **reads**: how many reads the value was read from. Below 20 nothing is
  reported as the value of the sample.
- **support**: the share of those reads carrying the value reported.
- **anchor support**: the share of the same reads agreeing on the
  residues where the anchor matched. The anchor is conserved, so this is
  how much the reads of this sample agree where nothing is difficult:
  the floor below which a value carried by few reads cannot be told from
  the errors of the reads.

Beside them the residues the reads carry where the anchor matched are
reported, with those of the reference. One residue of difference is a
strain difference and means nothing.

A read which shares no residue at all with the anchor of the reference
gives no value. The anchor is conserved, so a read which matched it
carries nearly all of it; a read which shares none of it matched bases
which translate to something else, because the read has an insertion or
a deletion there and the match is read a base out of frame, or because
the match is not of the anchor. Whatever lies beside such a match is
not the position the marker names. The residues that read carries are
still counted, so the disagreement stays visible in the report.

The test is made on each read, not on the sample: a sample can hold
both kinds of read at once, and the ones which carry the anchor still
give their value.

| verdict | meaning |
|---|---|
| `firm` | at least 60% of the reads agree on the value, and they agree at the anchor |
| `confirm` | the value is the majority, but too few reads agree for it to stand alone |
| `unresolved` | fewer than 35% of the reads agree on any value: the most frequent one is not the value of the sample, and none is reported |
| `too few reads` | fewer than 20 reads |

A value carried by at least a tenth of the reads and by at least 3
reads counts as one the sample holds. More than one such value is a mixture,
and is reported as such, with the share of each: this is where reading
the reads differs from reading a consensus.

The verdict is about the value, not about the virus. A `firm` verdict
means the reads agree; it does not mean the sample contains one virus.

### The markers read over a long window

`PB1-F2` and `NA-stalk` are read over far more than one codon: the
first needs the position of the first stop codon in a window of 92
codons, the second the distance between two anchors several hundred
bases apart. One indel anywhere in such a window changes what is read,
so on Nanopore reads they are often `unresolved` or read from too few
reads, and the distribution of their values is spread over neighbouring
lengths: on the samples measured below, `PB1-F2` is unresolved in 73 of
104 and `NA-stalk` has too few reads in 72. They are reported that way
rather than left out, because on Illumina reads, whose indel rate is
far lower, they are not in that state. `NS1-tail` is read over 16
codons but does not suffer from this, because its value is read
backwards from the stop codon, and is `firm` in 91 of the 104.

## Where the numbers come from

### The anchors, in the references flucut ships

Every anchor is checked against the whole built-in reference set when a
run starts, not only against the test data. In the influenza A set (136
records) the anchors are built as follows:

| marker | references of its segment | anchors built |
|---|---|---|
| `PB2-627`, `PB2-701` | 6 PB2 | 6, 6 |
| `PA-38` | 6 PA | 6 |
| `M2-31` | 6 M | 6 |
| `NS1-tail` | 6 NS | 6 |
| `PB1-F2` | 6 PB1 | 6 |
| `NA-275` | 36 NA | 34 |
| `NA-stalk` | 36 NA | 36 |
| `HA-226-228` | 64 HA | 56 |

The eight HA references without the receptor-binding anchor are the H13
and H16 ones, whose 220-loop lacks the conserved tryptophan; the two NA
references without the oseltamivir anchor are N7 ones, in which the
bases of the anchor are also found, within four errors, on the other
strand elsewhere in the segment. In neither case does an anchor land
somewhere else: the marker is simply not reported for a sample of that
reference.

In the influenza B set (32 records) **no anchor of the table is found at
all**, so an influenza B sample reports no marker rather than a wrong
one.

### The values, against the consensus of the same reads

The six markers read as single residues were compared with the IRMA
consensus of the same reads, on 104 Nanopore samples of the
multisegment RT-PCR (the set used for the cleavage site: 101 H5N1, 2
H5N2 and 1 H9N2, 94 R10 and 10 R9 samples, chosen for a high share of
defective genomes, so a hard case). The consensus was read at the same
anchored position, so the two readings are of the same residue. A
sample is compared where the reads gave a value and the anchor peptide
is in the consensus.

| marker | agree | differ |
|---|---|---|
| `PB2-627` | 87 | 1 |
| `PB2-701` | 88 | 0 |
| `PA-38` | 103 | 0 |
| `NA-275` | 99 | 0 |
| `HA-226-228` | 100 | 0 |
| `M2-31` | 97 | 0 |
| **total** | **574** | **1** |

The one difference is a sample whose reads carry both residues at PB2
627, E in 76% and K in 22% of 85 reads, where the consensus has to pick
one and picks the other.

### What the reads said

| marker | samples with a value | the values | `firm` | `confirm` | `unresolved` | too few reads |
|---|---|---|---|---|---|---|
| `PB2-627` | 97 | E 50, K 47 | 74 | 23 | – | 7 |
| `PB2-701` | 104 | D 104 | 93 | 11 | – | – |
| `PA-38` | 103 | I 103 | 98 | 5 | – | 1 |
| `NA-275` | 99 | H 99 | 93 | 6 | – | 5 |
| `HA-226-228` | 100 | QSG 77, QRG 22, LMG 1 | 98 | 2 | – | 4 |
| `M2-31` | 104 | S 93, N 11 | 102 | 2 | – | – |
| `NS1-tail` | 104 | ESEV* 98, EPEV* 2, other 4 | 91 | 13 | – | – |
| `NA-stalk` | 11 | 429 (9), 431 (2) | – | 11 | 21 | 72 |
| `PB1-F2` | 24 | 87 (9), 90 (6), 8 (4), 57 (4), 11 (1) | 7 | 17 | 73 | 7 |

The residues are those of avian H5N1: no oseltamivir and no baloxavir
residue, the avian receptor-binding pair in all but one sample, and the
avian PDZ-binding motif ESEV. PB2 627 is lysine in about half of them,
which the consensus of the same reads confirms in 87 of the 88 samples
where both could be read.

### What a consensus does not show

In 43 cases a residue is carried by at least a tenth of the reads and
is not the one in the consensus of those same reads: 30 at PB2 627, 6
at the receptor-binding site, 5 at M2 31 and 2 at PB2 701. The clearest
are at M2 31, where the reads are many: asparagine in 41% of 55,723
reads in one sample, in 35% of 41,774 in another, in 17% of 7,130 and
15% of 12,865 in two more, with serine in the consensus of each. Those
samples hold two M segments, one of which resists the adamantanes;
reading the position off the consensus alone shows only the other.

### How much the reads agree where nothing is difficult

The share of the reads agreeing on the residues where the anchor
matched, over all markers of all samples with at least 20 such reads:

| chemistry | markers | median | 10th percentile |
|---|---|---|---|
| R10 | 833 | 0.864 | 0.727 |
| R9 | 85 | 0.615 | 0.384 |

The two chemistries separate by themselves, as the control windows of
the cleavage site do, so the floor below which a minority is not told
from the errors of the reads is read from the reads in hand and not
from a rule naming the chemistry. It is also why a value below a tenth
of the reads is not counted as one the sample holds: on R9 reads the
disagreement at the anchor alone is of that size.

### How many reads reach a marker

The anchor is found in a median of 73% to 85% of the reads searched for
the markers near a segment terminus (`PB2-701`, `PA-38`, `M2-31`,
`NS1-tail`, `PB1-F2` and the first anchor of `NA-stalk`), and in 0.8%,
4.4% and 19% for `PB2-627`, `HA-226-228` and `NA-275`, which lie in the
middle of their segment. These samples were chosen for a high share of
defective genomes, which keep the two ends of a segment and lose the
middle, and the markers in the middle are also the ones whose reads an
internal priming site takes away. The read counts reported say so.

`NA-stalk` suffers most, because it needs both of its anchors, some 430
bases apart, in one read: 72 of the 104 samples have too few such reads
for any value.

### What a read out of frame at the anchor removes

Refusing a read which shares no residue with the anchor changed one
reported value in the 104 samples: at PB2 701 in one sample, from K to
D. There 178 of its 204 reads carried `TEDTDQH*` where the anchor
matched and the reference has `KRYGPALS` — the same bases read one base
out of frame — and the 26 reads which carried the anchor gave D, as
every other sample and every reference do. One further verdict went
from `confirm` to `firm`, and four of the minority residues above
disappeared.

It costs few reads: the median number of reads giving a marker falls by
at most 2% (`PB2-701` 1.8%, `NS1-tail` 2.0%, the others less). The
exception is a strain whose own anchor differs: in two samples the NS1
of the virus shares no residue with the anchor of any built-in
reference, so 96% of their reads are refused at `NS1-tail`. Both still
report ESEV* as `firm`, from 2,402 and 4,511 reads instead of 75,485
and 103,119.

The margin is wide enough for a strain difference: between any two
built-in references of the same marker the fewest shared of these eight
residues is two (the receptor-binding anchor, between the most distant
subtypes).

### What the kept part of a read removes

Reading a marker only from the bases flucut keeps changed the value
reported in 5 of the 746 sample-marker pairs with a value, and every
one of the five moved to the residue the consensus of the same reads
carries: four at PB2 627 from K to E, one at the receptor-binding site
from SKS to QSG. The disagreements with the consensus fell from 6 to 1.
Measured by running the same 104 samples on a build which reads the
markers from the whole read instead.

It costs reads only where the artefact was: the median number of reads
giving a marker is unchanged for six of the nine, 98% of what it was
for `PB1-F2` and 92% for `HA-226-228` and for `PB2-627`, whose anchor
sits on an internal priming site of the Uni12 primer (flucut reports
that site, at PB2 1913 to 1918, in its own table of internal priming).

## What it does not do

- It does not call a clade, a genotype or a lineage. A marker is a
  residue beside a conserved anchor, readable from one read; a clade is
  a placement on a tree, which needs a near-complete haemagglutinin and
  a nomenclature that is revised. That work belongs to a classifier run
  on the consensus afterwards — Nextclade with an avian A(H5) dataset,
  or LABEL, which IRMA ships — with the dataset tag recorded beside the
  result. The manual says the same in the user's terms
  (`docs/manual/markers.rst`).
- It does not predict a phenotype. A residue associated with resistance
  in the literature is reported as that residue, not as resistance; a
  susceptibility result needs a susceptibility assay.
- It does not replace a consensus. A marker is a handful of positions
  read early; the rest of the segment still has to be assembled.
- It says nothing about positions outside the table.

## Where it is reported

| output | what it holds |
|---|---|
| `flucut_summary.json` | the `markers` array of each group: every field of every marker, with the values and their shares |
| `summary.tsv` | one column per marker, with the value only where the reads support one, and `markers_attention` naming the markers which deserve a look |
| the HTML report | the **Markers** card of each sample, with the values, the mixtures and what each marker is; in the table of a segment across the samples, a column with the values of its markers |
| the dashboard | the same card and column, for the sample on the screen |
| the flags of a sample | one remark naming the markers whose residue deserves a look, with the share of the reads carrying it |

A residue which deserves a look is one the table marks as worth
attention: 627K, 701N, 38T, 275Y, 226L or 228S, and a stalk shorter than
the one of the reference by at least nine bases. The amantadine residue
M2 31N is reported but raises no remark: nearly every virus now
circulating carries it, so a remark on it would be a remark on every
sample.

## In the code

| file | what it holds |
|---|---|
| `src/seq/search/anchor-site.h`, `.c` | an anchor: building it from a reference, finding it in a read, reading the bases beside it. Not specific to the markers |
| `metadata/markers.yaml` | **the table itself**: one record per marker, with its anchor, its positions, their numbering convention and its meaning. A marker is added or changed here, not in the sources |
| `scripts/markers-table.awk` | turns that file into `obj/markers-table.h` when the program is built (awk only, so that building flucut still needs nothing but a C compiler and `make`) |
| `src/seq/markers/markers.h`, `.c` | the patterns, and the anchors of each marker in each reference |
| `src/seq/markers/marker-read.h`, `.c` | reading the markers out of one read |
| `src/analysis/marker-counts.h`, `.c` | the values counted over the reads of a sample, and what they say |
| `src/core/trimmer.c` | which reads are read, per platform, and the junction of a read with an internal deletion |
| `src/analysis/counter.c` | the values of a group of input files, merged over the threads and saved for `--resume` |
| `src/reports/json/json-markers.c`, `src/reports/html/html-report-markers.c`, `src/io/results.c` | the outputs |
