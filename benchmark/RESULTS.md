# Results

## 2026-10-30: low coverage (simulated reads)

What was run: `benchmark/run-coverage.sh` with the genome of
`make-genome.sh`: 10, 20 and 40 reads per segment, seeds 1 to 3 (other
parameters as in `simulate.py`), the tools of `run-tools.sh` with the
same commands (flucut at commit 3bfb35c), IRMA 1.4.0 FLU-minion. The
consensus compared is the majority consensus (`A_<gene>.fasta`,
`compare_consensus.py --plain`): the amended consensus of FLU-minion is
all N at these coverages, as it replaces with N the bases supported by
fewer than 50 reads (`MIN_CONS_SUPPORT=50`), whatever the trimming.

Differences from the truth (substitutions, insertions and deletions,
summed over the 8 segments) and the codon of the I348K site:

| reads per segment | seed | untrimmed | flucut | cutadapt | IRMA-core w70 |
|---|---|---|---|---|---|
| 10 | 1 | 4, **AAA** | 2, ATA | 2, ATA | 2, ATA |
| 10 | 2 | 36, **AAA** | 34, ATA | 34, ATA | 34, ATA |
| 10 | 3 | 0, ATA | 0, ATA | 2, ATA | 0, ATA |
| 20 | 1 | 0, ATA | 0, ATA | 0, ATA | 0, ATA |
| 20 | 2 | 1, **AAA** | 0, ATA | 0, ATA | 0, ATA |
| 20 | 3 | 1, ATA | 1, ATA | 1, ATA | 1, ATA |
| 40 | 1 | 1, **AAA** | 0, ATA | 0, ATA | 0, ATA |
| 40 | 2 | 0, ATA | 0, ATA | 0, ATA | 0, ATA |
| 40 | 3 | 1, **AAA** | 0, ATA | 0, ATA | 0, ATA |

All the runs assembled the 8 segments. The untrimmed reads give the
artefact in 5 of the 9 runs; every trimming tool removes it. The errors
which remain are the same for the three tools, from the low coverage
itself, except one run (10 reads, seed 3) in which cutadapt has 2 more:
on these reads, the template cut by cutadapt does not change the
consensus measurably, even at 10 reads per segment.

## 2026-10-30: flucut with truncated primers at the read ends

The change: at a read end, a primer is also searched without up to half
of its bases on the outer side (with one error less than its length
allows), and the outermost match near the end is the primer of that end
(before: the match with the fewest errors). Same reads and scoring as
below; flucut before (4e819ae) and after the change, default options:

| flucut | written | template removed | primer kept | flank kept | internal site kept |
|---|---|---|---|---|---|
| before | 3177 | 1280 | 3403 | 10247 | 31 |
| after | 3176 | 882 | 1852 | 7557 | 23 |
| (cutadapt, below) | 3177 | 102408 | 3037 | 8784 | 25 |

The read not written any more is a primer dimer (no template). A first
version allowed the truncated primers as many errors as their length
allows: it removed more template (1910 bases: 8-base truncated forms
with an error matched by chance inside the template of reads without
the primer), hence one error less. The time of the run went from 0.7 s
to 1.3 s (one thread).

## 2026-10-30: tools compared (simulated reads)

What was run: `benchmark/run-tools.sh` on the 3200 simulated reads of
the IRMA checks (`work/sim`), flucut at commit 4e819ae (default
options), cutadapt 5.2, IRMA-core 0.11.0 (trimmer), IRMA 1.4.0
(FLU-minion, as it is and with `SKIP_E=0`). The commands are in the
script:

- cutadapt: `-e 0.2 --revcomp -n 2 -m 1 -g AGCRAAAGCAGG
  -a CCTTGTTTCTACTGATCACGCGT` (both primers, not anchored, as the
  reads have flanking bases; both orientations);
- IRMA-core trimmer: the two forms of Uni12 and the tailed Uni13 as
  primers, `--p-kmer-length 12` (the length of Uni12; MIRA-NF uses 14
  for custom primers, longer than Uni12), `--p-fuzzy`, windows
  `--p-restrict` 40 (the MIRA-NF default for custom primers) and 70.
  MIRA-NF itself does no primer trimming for influenza Nanopore reads
  ("Flu-ONT experiment type does not need primers").

Reads (bases by their true origin, see `score_reads.py`):

| tool | written | lost | template removed | primer kept | flank kept | internal site kept |
|---|---|---|---|---|---|---|
| untrimmed | 3200 | 0 | 0 | 102554 | 127032 | 2698 |
| flucut | 3177 | 0 | 1280 | 3403 | 10247 | 31 |
| cutadapt | 3177 | 0 | 102408 | 3037 | 8784 | 25 |
| IRMA-core, window 40 | 3195 | 0 | 77 | 43474 | 75068 | 1105 |
| IRMA-core, window 70 | 3200 | 0 | 130 | 13152 | 22890 | 660 |

IRMA consensus, summed over the 8 segments (`_elong`: with the
reference elongation); the codon of the I348K site in the majority
consensus of HA:

| run | missing 5' | missing 3' | beyond 5' | beyond 3' | errors | ambiguous | codon |
|---|---|---|---|---|---|---|---|
| untrimmed | 219 | 276 | 0 | 0 | 0 | 2 | **AAA** |
| untrimmed_elong | 0 | 1 | 140 | 208 | 11 | 10 | **AAA** |
| flucut | 219 | 276 | 0 | 0 | 0 | 0 | ATA |
| flucut_elong | 96 | 104 | 0 | 1 | 0 | 0 | ATA |
| cutadapt | 219 | 276 | 0 | 0 | 0 | 0 | ATA |
| cutadapt_elong | 96 | 104 | 0 | 1 | 0 | 0 | ATA |
| IRMA-core w40 | 219 | 276 | 0 | 0 | 0 | 1 | ATA |
| IRMA-core w40_elong | 0 | 0 | 128 | 218 | 18 | 10 | ATA |
| IRMA-core w70 | 219 | 276 | 0 | 0 | 0 | 0 | ATA |
| IRMA-core w70_elong | 84 | 102 | 0 | 1 | 11 | 4 | ATA |

(The 1 base beyond the 3' end with the elongation is an N appended by
IRMA to segment 3.)

What this shows, on these reads:

- **The consensus**: cutadapt is as good as flucut (the same consensus,
  without errors, 12 and 13 bases short at the ends with the
  elongation). The IRMA-core trimmer leaves primers and flanks in many
  reads (its 12-mers need an exact or one-mismatch match, which
  Nanopore errors often break): with the elongation its consensus has
  errors, and with window 40 it extends into flanks. Every trimming
  tool removes the I348K artefact: the primer bases left at the site
  (1105 of 2698 at most) stay a minority.
- **The reads**: cutadapt removes 80 times as many template bases as
  flucut (102408 against 1280): its non-anchored 5' primer cuts
  wherever Uni12 matches with 2 errors, also deep inside a read; 166
  reads lose more than 50 bases each, among them 21 complete HA reads
  cut at the I348K site (about 1100 bases each), and reads of PA, PB1
  and PB2. With 400 reads per segment this does not change the
  consensus; it reduces the coverage, which matters at low coverage
  and for the frequencies of the variants (not measured yet).
- **flucut's own errors**: the 1280 template bases are mostly from 15
  internally primed HA reads (about 68 bases each) which begin inside
  the primer: the partial primer was not found, and a Uni12-like
  sequence 55 bases after the site (HA 1145-1157, `AGCAACGAGCAGG`,
  edit distance 2) was taken as the primer, as it is within the window
  of 100 bases. At the internal site, flucut leaves 31 primer bases,
  cutadapt 25.

### Barbell

Not run: its custom mode does not fit a single primer per end without
barcodes. With one query per end it stops ("For now we only support
'groups' if you only have one query, just add a 'random' second query");
with the two forms of Uni12 as the left queries and a decoy second
right query differing in one base (its documented workaround), it set
the flank error cutoffs to 0 (left) and 4 (right) and annotated 66
matches in the 3200 reads (none on the 12 exact reads of
`testdata/reads.fastq`); with `--flank-max-errors` 1 to 4, at most 49
reads with the right primer, while the left flank (11 bases) matched
all over the reads at 4 errors (7038 matches). Barbell finds the
flanks around a barcode (adapter and primer) and is meant for barcoded
reads: the fair comparison needs reads simulated with the adapters and
barcodes of a kit (from the kit definitions of Barbell), with
`adapter + barcode + primer` queries, and `--use-extended` (described
as applying to kits). Barbell 0.3.3, commit 057a902.

## 2026-10-30: read-level scores of flucut (simulated reads)

`benchmark/score_reads.py` on the reads of the IRMA checks (3200 reads;
`work/sim`), flucut at commit 45ceab0, default options, and with
`--trim terminus -r` the genome. Bases by their true origin (see
`score_reads.py`):

| reads | written | lost | template removed | primer kept | flank kept | internal site kept |
|---|---|---|---|---|---|---|
| untrimmed | 3200 | 0 | 0 | 102554 | 127032 | 2698 |
| flucut | 3177 | 0 | 1280 | 3403 | 10247 | 31 |
| flucut --trim terminus | 3177 | 0 | 1035 | 73185 | 10681 | 31 |

The 23 reads not written are the primer dimers (no template). The
flanking bases kept are mostly at read ends without a primer (reads
broken in the middle, reads with one primer), where nothing can be
trimmed by position. The 31 bases of the internal site kept come from
internally primed reads whose Uni12 primer was not found (more errors
than allowed). With `--trim terminus`, the conserved termini are kept
by design (primer kept).

## 2026-10-30: length check (simulated reads)

flucut with `-r` the genome, on the reads of the IRMA checks: for the
HA (KF918483.1, 1779 bases), 209 of 340 complete amplicons are shorter
than 90% of the reference, most at 670-679 bases (the length from the
internal site at 1078 to the end of the segment, between the primers:
1779 - 1077 - 12 - 13 = 677); in the other segments, at most 1 (PB2:
1 of 338, at 730-739). Without references, the log lists the peaks
2310-2319, 860-869, 1000-1009, 1380-1389, 2200-2209, 1540-1549,
670-679 and 1750-1759: the eight segments, and the internal amplicons
of HA, more frequent than the complete HA (152 against 80 reads in
those bins), as the simulation starts 60% of the HA reads at the site.

## 2026-10-30: speed of the classification against references

What was run: flucut on the 3200 simulated reads of the IRMA checks
below (`work/sim/reads.fastq.gz`), one thread (`-j 1`), with no
references, with the 8 segments of the genome (`-r genome.fasta`), and
with 80 references (the genome 10 times, renamed); with the k-mer index
(default) and with `--exhaustive-search`. Wall-clock time of the whole
run, two runs each, on the 4-core machine of the benchmark.

| references | no `-r` | k-mer index | exhaustive |
|---|---|---|---|
| 8 | 0.68-0.72 s | 0.76-0.79 s | 2.0-2.2 s |
| 80 | | 1.1-1.2 s | 14.0-16.4 s |

The outputs with the index and exhaustive are identical (compared
byte by byte, decompressed), with 8 and with 80 references. Without the
exhaustive search for the read ends without a candidate, 47 read ends
(of about 5800) were not located by the index (reads with an error in
each 10-mer of the probe); with it, none.

## 2026-10-30: IRMA checks (simulated reads)

What was run: `benchmark/run-irma-checks.sh` (see `README.md`), on the
genome of `make-genome.sh` (downloaded 2026-10-30), with the defaults of
`simulate.py` (seed 1, 400 reads per segment, 3200 reads, 239 of them
started at the internal site of HA).

Versions: flucut at commit 4eec29e (with the benchmark scripts of this
commit; the program is unchanged); IRMA 1.4.0 (commit 1311c81),
IRMA-core 0.11.0, IRMA-viz 0.2.2, LABEL 0.7.2; module FLU-minion.

Runs: `raw` (untrimmed), `flucut` (default trimming), `notags` (the same
with `--no-tags`), `terminus` (`--trim terminus -r` the genome); each
with FLU-minion as it is and with `SKIP_E=0` (`_elong`).

### The I348K site

The codon of the site in the plain (majority) consensus of HA:

| run | codon |
|---|---|
| raw, raw_elong | **AAA** (Lys: the artefact) |
| flucut, notags, terminus (both modes) | ATA (Ile: the truth) |

Untrimmed, at the middle base of the codon 61% of the reads carry the A
of the primer and 37% the T of the template (IRMA allele table of the
raw run); the majority consensus calls AAA, the amended consensus an
ambiguity (W). With flucut the site is the truth in all runs. With 0.6
of the HA reads internally primed (a simulation parameter, not a
measured share), the untrimmed reads reproduce the published artefact.

### Tags in the read titles

`flucut` and `notags` give the same IRMA read counts and the same
amended consensus in all 8 segments, with both configurations: IRMA
accepts the tags. `--no-tags` stays an option, not the default.

### The segment ends

Bases of the true segment missing (or beyond it) at each end of the
amended consensus, and errors:

| run | 5' end | 3' end | errors |
|---|---|---|---|
| any input, FLU-minion as it is | 20-45 missing | 20-58 missing | none (raw: 2 ambiguous at the HA site) |
| raw_elong | 0-28 **beyond** the truth | 13-32 beyond | 1-6 insertions per segment, ambiguous bases |
| flucut_elong | 12 missing | 13 missing | none (segment 3: one N appended at the 3' end) |
| terminus_elong | complete | complete | 1 ambiguous base per segment; segment 2: 1 substitution; segments 3 and 4: one N appended |

- With FLU-minion as it is, IRMA does not extend the consensus beyond
  its references: the same 20-58 bases are missing at the ends whatever
  the trimming, so keeping the termini changes nothing.
- With the elongation, the consensus of the flucut reads ends exactly
  where the primers began (12 and 13 bases short), without errors; the
  untrimmed reads extend it into flanks and primer tails, with
  insertions.
- Keeping the termini (`--trim terminus`) gives complete ends, but the
  bases there are the primer's: the 4th base is an ambiguity (R) in all
  segments, from the degenerate Uni12 primer (the even mix of A and G
  is a simulation assumption), and in segment 2 the 3' terminus of PR8,
  `CCTTGTTCCTACT`, is overwritten by the primer, `CCTTGTTTCTACT` (real
  sequences, no assumption). This supports removing every primer by
  default.

### Not yet measured

The other tools of `TOOLS.md`, the read-level measures (primer bases
left, template bases removed, reads lost: the truth is in
`truth.tsv.gz`), other shares of internal priming, other error rates,
and real data.
