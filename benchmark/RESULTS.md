# Results

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
