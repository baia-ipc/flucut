# Benchmark

What the benchmark measures, on which data, and how to run it. The
results, with the versions and the commands, are in `RESULTS.md`; the
tools compared in `TOOLS.md`.

The measure which decides is the IRMA consensus: whether a variant of
the consensus comes from a primer, at the segment ends or at an
internal priming site. The reads lost and the bases left or removed at
the read ends come next; speed and memory last.

## Data: simulated reads

For now only simulated reads are used, so that the truth is known for
every base.

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
becomes the A of `AAA` (Lys). This is the I348K artefact of
`docs/internal-priming.md`, found in the public sequences. In
KF918483.1 it is the only position with as few mismatches (2), at
1078-1089 (1-based).

### The reads

`simulate.py` simulates amplicons of the default scheme of flucut (the
degenerate untailed Uni12 primer and the tailed Uni13 primer), 400 reads
per segment by default. A share of the HA reads (default 0.6) starts at
the internal site with the bases of the primer. Each read gets an
orientation, random flanking bases (0-40 at each end, standing for
adapter and barcode remnants: random, not a real adapter), sometimes a
start or an end inside a primer (10% each) or a break in the middle
(5%), and substitutions, insertions and deletions at a rate drawn for
each read (mean 4%). 1% of the full amplicons are primer dimers. The R
of Uni12 is A or G at random in each read. All of these are modelling
choices, not measured properties of a run. The origin of each base of
each read is written to `truth.tsv.gz`.

## IRMA

IRMA 1.4.0 with IRMA-core 0.11.0, IRMA-viz 0.2.2 and LABEL 0.7.2,
installed from the sources by `setup-irma.sh` (the release archives were
not reachable from the environment). The module is `FLU-minion` (the
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
