# Benchmark

What the benchmark measures, on which data, and how to run it. The
results, with the versions and the commands, are in `RESULTS.md`; the
tools compared in `TOOLS.md`.

The measure which decides is the IRMA consensus: whether a variant of
the consensus comes from a primer, at the segment ends or at an
internal priming site. The reads lost and the bases left or removed at
the read ends come next; speed and memory last.

"The choice of the reference of each segment" below is a study of its
own: it measures how flucut picks the reference of each segment, and
whether a mixed infection can be seen in the votes. It has its own
genomes, scripts and measure (the read ends classified correctly).

## Data: simulated reads

Most measures use simulated reads, so that the truth is known for
every base. Four public sequencing runs are measured as well ("Public
runs" below).

### The genome

`make-genome.sh` downloads it from NCBI:

- segments 1-3 and 5-8: A/Puerto Rico/8/1934 (H1N1), RefSeq
  NC_002016.1 to NC_002023.1 (except NC_002017.1, its HA): complete,
  with both conserved termini;
- segment 4: the HA of A/chicken/Cambodia/X0815302/2013 (H5N1),
  KF918483.1, complete with both termini.

The combination (PR8 internal genes with an H5 HA) is constructed for
the simulation; it is not a sequenced isolate.

Why this HA: in the Cambodian H5N1 HA sequences of GenBank (30 checked),
the bases after the motif GLFGAIAG read `AGCTATAGCAGG`, which differs
from the Uni12 primer `AGCRAAAGCAGG` at two positions: the R is on the
third base of an Ala codon (silent), and the T of the Ile codon `ATA`
becomes the A of `AAA` (Lys). This is the I351K artefact of
`docs/internal-priming.md`, found in the public sequences (residue 6 of
HA2, at the start of the fusion peptide; 351 is counted on the HA0
precursor from the initiator methionine). In
KF918483.1 it is the only position with as few mismatches (2), at
1078-1089 (1-based).

### The reads

`simulate.py` simulates amplicons of the default scheme of flucut for
the type of the genome (`--virus`): for influenza A the degenerate
untailed Uni12 primer and the tailed Uni13 primer, the same in every
segment; for influenza B the pair of the IBV-GA2 cocktail designed for
each segment. 400 reads per segment by default. A share of the HA reads (default 0.6) starts at
the internal site with the bases of the primer. Each read gets an
orientation, random flanking bases (0-40 at each end, standing for
adapter and barcode remnants: random, not a real adapter), sometimes a
start or an end inside a primer (10% each) or a break in the middle
(5%), and substitutions, insertions and deletions at a rate drawn for
each read (mean 4%). 1% of the full amplicons are primer dimers. The R
of Uni12 is A or G at random in each read. All of these are modelling
choices, not measured properties of a run. The origin of each base of
each read is written to `truth.tsv.gz`.

With `--wrong-pair-share` (influenza B), a share of the amplicons began
at the terminus of their segment but were primed by the pair of other
segments: their primer bases are not the sequence of that segment, and
flucut has to remove all of them. The evaluation of influenza B
(`../influenza_b_evaluation/`) measures that.

With `--joined-share`, a share of the reads is two amplicons joined end
to end (the second one of any segment, each in its own orientation,
0-20 other bases between them): a tool which splits reads is scored on
them, by the reads split wrongly and those not split
(`score_reads.py`). The default has none.

## IRMA

IRMA 1.4.0 with IRMA-core 0.11.0, IRMA-viz 0.2.2 and LABEL 0.7.2: an
installed release, or one built from the sources by `setup-irma.sh`.
The module is `FLU-minion` (the
configuration of IRMA for Nanopore reads), as it is, and with the
reference elongation on (`SKIP_E=0`), which extends the consensus into
the untranslated regions of the segments.

## Running

    benchmark/make-genome.sh genome
    benchmark/setup-irma.sh tools
    benchmark/run-irma-checks.sh genome/genome.fasta work bin/flucut \
        tools/irma/IRMA

`run-irma-checks.sh` simulates the reads, trims them, runs IRMA and
compares each consensus with the truth (`compare_consensus.py`).
`score_reads.py` scores the trimmed reads of any tool against the truth
of each base:

    python3 -I benchmark/score_reads.py work/sim work/flucut/reads.fastq.gz \
        flucut

It counts the reads written and lost, the template bases removed, the
primer bases kept (in all, and apart those of the primer tails and
those of the conserved termini, which flucut keeps on purpose at the
terminal ends when it has references), the flanking bases kept, the
primer bases kept at the internal site, and the reads split: flucut
writes a read which holds several amplicons as several reads, so it can
write more reads than it was given.

flucut is run in its configurations side by side, as its result
depends on the references: with its default options (the reads are
classified against the built-in references and the conserved termini
kept at the terminal ends), with `-r off` (no references: every primer
base removed; the primers of influenza A and B are then both searched)
and with `-r off --virus a` (the primers of influenza A alone).

The other scripts, each with its usage in its first lines; the
programs are given in the environment (`FLUCUT`, `IRMA`, `IRMA_CORE`,
`BBDUK`, `PORECHOP`, `DORADO`, `BARBELL`, `SAMTOOLS`, `MINIMAP2`;
cutadapt and a `python3` with edlib in the `PATH`), and a tool whose
program is not given is left out:

| script | what it does |
|---|---|
| `run-tools.sh <simulation> <outdir>` | every tool on one set of simulated reads: read-level scores, IRMA consensus with and without the elongation, the codon of the site |
| `run-coverage.sh <genome> <workdir> <flucut> <IRMA> <irma-core>` | 10, 20 and 40 reads per segment, three seeds each: flucut, cutadapt and the IRMA-core trimmer, IRMA majority consensus |
| `run-flucut-options.sh <simulation> <outdir>` | flucut alone, read-level scores by its options: references, adapter search, errors allowed in a primer |
| `run-speed.sh <simulation> <outdir> [<copies>]` | time and memory of every tool (`measure.py`) on the reads repeated 10 times |

The reads of the other sets are simulated with the options of
`simulate.py`: `--barcodes 4` (native barcoding adapters outside the
primers), `--error-rate 0.01` to `0.12`, `--joined-share 0.05`.

`measure.py` reports the wall time, the CPU time and the peak memory
of a command; with `MEASURE_CPUS=<list>` in the environment the command
is confined to these processors (e.g. 4 cores of a larger machine).

## Public runs

`run-real-nanopore.sh` and `run-real-illumina.sh` take a run accession
of SRA, the accessions of the eight segments of its isolate in GenBank,
an output directory and, optionally, a directory for the downloads:

    benchmark/run-real-nanopore.sh SRR36760543 \
        PZ239631.1,PZ239632.1,PZ239633.1,PZ239634.1,PZ239635.1,PZ239636.1,PZ239637.1,PZ239638.1 \
        work/SRR36760543 reads/SRR36760543

They download the reads (ENA) and the segments (NCBI), trim the reads
with each tool (timed), map them to the segments of the isolate and
count the codon of the internal priming site by where the reads start
(`site_reads.py`), and assemble the reads of each tool with IRMA; the
consensus sequences are compared with that of the untrimmed reads.
flucut is run with its default options, with `-r off` and with the
segments of the isolate as references; for the Nanopore runs also
with the file given in `REFERENCES` (complete segments, with the
untranslated regions, which most GenBank records lack). `IRMA_JOBS`
sets how many assemblies run at a time.

`read_origins.py` writes, from the JSON summaries of flucut runs, the
table of the 5' read ends by origin (terminal, internal, not located)
and the largest internal priming sites:

    python3 -I benchmark/read_origins.py \
        SRR36760543=work/SRR36760543/flucut/flucut_summary.json

## The choice of the reference of each segment

flucut picks one reference per segment for each sample, by the votes of
probes taken from the first reads (`src/auto-settings.h`). Three
constants decide: how many votes a reference needs
(`AUTO_REFERENCE_MIN_VOTES`), the share of its segment's votes from
which its lead is clear (`AUTO_REFERENCE_MIN_SHARE`), and the share of
the leading reference's votes from which another one is close enough to
be followed (`AUTO_REFERENCE_CLOSE_SHARE`).

The votes alone cannot say which value is right: the vote share says
how similar the reads are to each reference, not which reference
classifies them better. This simulation study answers it by simulating
reads of genomes whose true sequence is known, and measuring what each
reference gets right. The results are in `RESULTS.md` ("the choice of
the reference of each segment").

### The genomes

`reference_choice_genomes.tsv` lists 17 complete genomes (the eight
segments, each with both conserved termini) of influenza A viruses
which are **not** in the built-in reference set, at graded distances
from it: two human seasonal viruses of other decades, a 1957 H2N2, two
South-East Asian H5N1 of poultry, a Chinese H7N9 of a human case (a
known reassortant: H7 and N9 with the internal genes of poultry H9N2),
poultry H9N2 and H10N8, clade 2.3.4.4 H5N8 and H5N6, a Eurasian and an
American wild-bird virus, equine and canine lineages, and three swine
viruses. Several have segments which sit between two references of the
set: the H7N9 and the swine viruses by their history, the wild-bird
viruses by their internal genes.

The table gives, for each genome, a label, the strain, the subtype, the
accessions of its eight segments and why it was chosen.
`fetch_genomes.py` downloads them from NCBI (the requests carry no
personal details, only a generic client name) and refuses a genome
whose segments do not carry both termini, as the simulation needs
complete segments:

    python3 -I benchmark/fetch_genomes.py \
        benchmark/reference_choice_genomes.tsv work/genomes

Complete termini are rare in GenBank: most records are coding
sequences. This is why the set is weighted towards older submissions,
and why no recent (2024-2025) cattle H5N1 genome is in it.

### How a sample is simulated

Each sample is simulated by `simulate.py` (see "Data: simulated reads"
above for the read model), with 400 reads per segment, the reads in
random order (`--shuffle`), as flucut takes its probes from the first
reads of a file. The error rate per read is drawn around a mean of 4, 8
or 12%: the range of the real Nanopore runs measured in `RESULTS.md`.

The truth of every base is written by the simulator, so for each read
end with a primer it is known whether the amplicon started at a
terminus of the segment or at the internal priming site, and at which
position of the simulated segment.

### What is measured

For a read end, flucut reports the origin (`ol`, `or`: terminal,
internal or unknown) and where its probe mapped (`ml`, `mr`:
reference, position, strand). `score_read_ends.py` counts a read end if
the simulated read has a primer there and flucut found it (what happens
before the classification does not depend on the reference), and counts
it correct if

- the origin is the true one, and
- the position is within 5 bases of the true one, on a reference of the
  right segment.

The true position is a position of the *simulated* segment; the
reference differs from it by substitutions, insertions and deletions.
It is carried over by a global alignment of the two (`pairwise.py`,
minimal edit distance), which also gives their identity. The score of a
segment is the share of its read ends which are correct.

The references passed to flucut may hold one reference of the segment
or several (flucut keeps them all where it cannot choose): with several,
a read end is correct if it is right on the reference it was mapped to.

### The design

`reference_choice.py` runs, for each of the 17 genomes and each of the
three error rates:

- a sample of 3200 reads (400 per segment), classified six times: once
  against the references of each of the five subtypes, and once against
  the whole set. This gives the accuracy of every candidate reference
  of every segment, and of keeping all of them;
- the same sample against the whole set, with the votes of every
  reference and the reference chosen taken from the log of flucut;
- the first 30, 60, 125, 250, 500 and 1000 reads of that sample, and
  two samples with few reads of some segments (6 to 50 reads of four
  segments), for the constants which decide when a segment has enough
  votes.

That is 51 samples, 2040 scored classifications and 3672 segment
decisions. The simulated reads are written below `<workdir>/samples`
and removed once scored; `accuracy.tsv` and `votes.tsv` remain.

    python3 benchmark/reference_choice.py work/genomes \
        ~/.cache/flucut/flu-a-references-v1.fasta bin/flucut work/choice 8

(the last argument is the number of samples simulated at a time.)

`reference_rule.py` is the rule of flucut (`auto_settings_decide_
references()` in `src/auto-settings.c`) with its three constants as
parameters, so that other values can be applied to the votes of a run
without rebuilding the program. It reproduces the choices of flucut in
all 3672 decisions of the study.

`reference_choice_tables.py` writes the tables of `RESULTS.md` from
`accuracy.tsv` and `votes.tsv`; the last argument picks one:

    python3 benchmark/reference_choice_tables.py work/choice \
        ~/.cache/flucut/flu-a-references-v1.fasta follow

| table | what it shows |
|---|---|
| `segments` | every simulated segment: the vote share, identity and accuracy of each reference, and the one chosen |
| `follow` | where following the common subtype changes the reference, that reference against the leading one, by the ratio of their votes |
| `rule` | what other values of the three constants choose, and the accuracy lost |
| `all` | the chosen reference against keeping all the references of the segment |
| `votes` | whether a leading reference with few votes is the one which leads with all the reads |
| `termini` | the bases each reference lacks at its ends, and the accuracy with it |

### Naming the subtype of a sample

flucut names the subtype of a sample after the references it chooses
for the HA and the NA, so a sample whose HA or NA is far from every
reference of its own subtype is named after another one.
`subtype_calls.py` measures how often that happens: it simulates one
sample per genome and error rate, runs flucut against the references,
and compares the subtype reported with the subtype of the genome. The
HA and the NA are counted separately, as they are chosen separately.

    python3 benchmark/subtype_calls.py work/genomes \
        benchmark/reference_choice_genomes.tsv \
        ~/.cache/flucut/flu-a-references-v3.fasta bin/flucut \
        work/subtypes 8 200

`subtype_panel_genomes.tsv` is a second table of genomes for the same
measure, which covers every HA subtype (H1 to H16) and every NA subtype
(N1 to N9): the HA of one subtype and the NA of one subtype with the
internal genes of A/Puerto Rico/8/1934, a combination put together for
the measure, as only the HA and the NA decide the subtype named. None
of its HA and NA records is in the built-in set. It is fetched and used
like the other one:

    python3 -I benchmark/fetch_genomes.py \
        benchmark/subtype_panel_genomes.tsv work/panel
    python3 benchmark/subtype_calls.py work/panel \
        benchmark/subtype_panel_genomes.tsv \
        ~/.cache/flucut/flu-a-references-v3.fasta bin/flucut \
        work/panel_subtypes 8 200

(the last two arguments are the samples simulated at a time and the
reads per segment). `subtype_calls.tsv` holds a line per genome and
error rate: the subtype reported, the references chosen for the HA and
the NA with their subtypes, and whether the H and the N are right.

### Mixed samples

A segment whose votes are split has two causes which the votes alone do
not separate: one strain lying between two references (the probes of
the same read vote for different references), or two strains in the
sample (the probes of a read agree, the reads disagree with each
other). flucut reports, for the HA and the NA, the reads of a second
subtype among the first reads of a sample, in the HTML report and the
JSON summary (`src/subtype-reads.h`), and the agreement of the probes
within a read for each segment and reference in the messages (`-v`,
`src/read-agreement.h`).

`mixed_samples.py` simulates, at the three error rates and with two
seeds each:

- a sample of each of the 17 genomes alone, which gives what is seen
  without a mix, including the genomes which lie between two
  references;
- mixes of two strains, with the minor one at 50, 10, 5, 2 and 1% of
  the reads: eight pairs of different subtypes (one of them two swine
  viruses which share their internal genes and the subtype of their
  NA) and, as a control, two H5N1 of the same clade.

Each strain is simulated by itself, with its own read names
(`--read-prefix`) and the number of reads of its share; `mix_reads.py`
joins the simulations in an order of its own for each sample. By
default a sample has 400 reads per segment, 3200 in all, so a minor
fraction of 1% is 4 reads per segment; flucut counts the probes of the
first 2000 reads. The reads per segment can be given, and the reads of
the NA: with a scarce NA, flucut counts the probes of up to 20,000
reads, ten times more reads of the HA.

    python3 benchmark/mixed_samples.py work/genomes \
        ~/.cache/flucut/flu-a-references-v3.fasta bin/flucut work/mixed 8
    # 32,000 reads, 40 of them of the NA
    python3 benchmark/mixed_samples.py work/genomes \
        ~/.cache/flucut/flu-a-references-v3.fasta bin/flucut \
        work/mixed_deep 8 4000 40
    python3 benchmark/mixed_samples_tables.py work/mixed one_subtype
    python3 benchmark/mixed_samples_tables.py work/mixed by_fraction

`mixed_samples_tables.py` writes, from `second_subtype.tsv` (what the
JSON summary of each sample says about the HA and the NA):
`one_subtype` (the samples with one subtype in the gene, and how many
are above the threshold of flucut), `remarked` (for each pair and minor
fraction, the runs above it with the right second subtype) and
`by_fraction` (the same, all pairs together); and from `agreement.tsv`
(the messages), for the eight segments by reference: `single` (the
statistic in the samples of one strain, whose highest value per segment
is the threshold), `detected` (the runs in which each mix is above it)
and `lowest` (the lowest minor fraction above it in every run).

### What the simulation assumes

The read model is the one of `simulate.py` (orientation, flanking
bases, reads starting or ending inside a primer, breaks, primer
dimers, the error rate drawn per read): modelling choices, not measured
properties of a run. Beyond them, this study assumes that

- the sample holds one virus (or two, in the mixes), with no
  contamination and no barcode crosstalk: a minor strain of about 1%
  cannot be told from crosstalk using the reads alone;
- every segment is amplified: the samples with few reads of some
  segments model uneven amplification, not a segment missing
  altogether;
- the true genome is one of the 17 chosen. They cover the distances to
  the built-in set which matter here, but they are not a sample of
  what a surveillance laboratory sees;
- the measure is the classification of the read ends. What the choice
  of reference does to the consensus of an assembler is not measured.
