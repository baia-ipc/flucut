# Results

## 2026-10-30: the inward extension of the primer alignment

What was run: flucut at commit 16d001d, as it is and without the
extension of the primer alignment inward at equal cost (the loop after
the choice of the end in `primer_search_end()`, removed for this test
only), scored with `score_reads.py` (template removed, primer kept,
flank kept, internal site kept).

| reads | as it is | without the extension |
|---|---|---|
| IRMA checks (4% errors) | 680, 1819, 7542, 23 | 472, 2136, 7546, 28 |
| barcoded | 267, 1019, 698, 44 | 92, 1311, 701, 55 |
| 12% errors | 5916, 3715, 7662, 101 | 5299, 4585, 7680, 125 |
| 1% errors | 607, 1073, 6457, 1 | 556, 1127, 6421, 5 |

Without the extension, fewer template bases are removed, and more
primer bases (and bases of the internal site) are kept, in every set:
the extension is kept.

## 2026-10-30: time and memory

What was run: `benchmark/run-speed.sh` (the commands of
`run-tools.sh`, timed by `benchmark/measure.py`: wall time, CPU time of
all the processes, peak resident memory of the largest process) on the
reads of the IRMA checks repeated 10 times (32000 reads, 53 million
bases, 43 MB gzipped), with 4 threads where the tool has an option for
it; flucut at commit f7553e7; a machine with 4 cores, otherwise idle.
Barbell on the barcoded reads (`simnb`) repeated 10 times, with flucut,
cutadapt and dorado on them as well.

| tool | wall s | CPU s | peak MB |
|---|---|---|---|
| flucut, 1 thread | 6.68 | 6.42 | 27 |
| flucut, 4 threads | 2.43 | 7.59 | 60 |
| flucut, 4 threads, 8 references (`-r`) | 2.48 | 7.61 | 78 |
| cutadapt, 4 cores | 3.71 | 12.61 | 41 |
| IRMA-core trimmer (1 thread) | 7.47 | 7.19 | 11 |
| BBDuk hdist 2, 4 threads (two passes) | 3.96 | 11.81 | 954 |
| Porechop_ABI tuned, 4 threads | 33.62 | 80.35 | 253 |
| dorado trim, 4 threads | 8.08 | 10.65 | 479 |
| minimap2 + samtools ampliconclip, 4 threads | 15.01 | 25.05 | 125 |

Barcoded reads:

| tool | wall s | CPU s | peak MB |
|---|---|---|---|
| flucut, 1 thread | 7.03 | 6.70 | 27 |
| flucut, 4 threads | 2.50 | 7.89 | 62 |
| cutadapt, 4 cores | 3.77 | 12.95 | 40 |
| dorado trim, 4 threads | 7.71 | 10.59 | 502 |
| Barbell, 4 threads (annotate, filter, trim) | 2.41 | 5.40 | 44 |

The gzip compression of the output (zlib, level 6) is a large part of
the time of flucut: with level 1 (rebuilt for this test only), 4.32 s
wall and 4.03 CPU s with 1 thread, 1.72 s wall with 4, for an output 32%
larger (54 instead of 41 MB); decompressing the input alone takes 0.8 s,
recompressing it with `gzip -6` 7.3 s. Barbell and cutadapt were given
the same gzipped input; Barbell writes uncompressed FASTQ.

## 2026-10-30: the errors in a primer, by error rate of the reads (-e auto)

What was run: `benchmark/simulate.py` with `--error-rate` 0.01, 0.02,
0.04, 0.08 and 0.12 (the mean of the error rate drawn for each read;
the genome of the IRMA checks, seed 1, 400 reads per segment, without
barcodes; 0.04 gives the reads of the other sections), then flucut at
commit dddfb6c with `-e` from 0.1 to 0.35 (`-j 4`), scored with
`score_reads.py`. With Uni12 (12 bases), 0.1 and 0.15 allow 1 error, 0.2
2, 0.25 and 0.3 3, 0.35 4; with mbtuni-13 (23 bases), 2, 3, 4, 5, 6, 8.

| mean error rate | -e | template removed | primer kept | flank kept | internal site kept |
|---|---|---|---|---|---|
| 0.01 | 0.1 | 56 | 2486 | 9017 | 85 |
| 0.01 | 0.15 | 57 | 1804 | 8239 | 85 |
| 0.01 | 0.2 | 607 | 1073 | 6457 | 1 |
| 0.01 | 0.25 | 1734 | 483 | 4009 | 1 |
| 0.01 | 0.3 | 1876 | 409 | 3576 | 1 |
| 0.01 | 0.35 | 5693 | 130 | 1070 | 13 |
| 0.02 | 0.1 | 165 | 4034 | 11744 | 114 |
| 0.02 | 0.15 | 172 | 2668 | 10114 | 114 |
| 0.02 | 0.2 | 842 | 1419 | 7451 | 11 |
| 0.02 | 0.25 | 2437 | 622 | 4458 | 5 |
| 0.02 | 0.3 | 2598 | 600 | 4136 | 5 |
| 0.02 | 0.35 | 7370 | 272 | 1306 | 5 |
| 0.04 | 0.1 | 196 | 9283 | 16898 | 227 |
| 0.04 | 0.15 | 224 | 5197 | 13283 | 227 |
| 0.04 | 0.2 | 680 | 1819 | 7542 | 23 |
| 0.04 | 0.25 | 2105 | 801 | 4824 | 3 |
| 0.04 | 0.3 | 2366 | 810 | 4601 | 3 |
| 0.04 | 0.35 | 7912 | 360 | 1400 | 3 |
| 0.08 | 0.1 | 438 | 26381 | 37011 | 548 |
| 0.08 | 0.15 | 531 | 16817 | 28099 | 548 |
| 0.08 | 0.2 | 1296 | 6540 | 14094 | 157 |
| 0.08 | 0.25 | 3609 | 2123 | 6519 | 64 |
| 0.08 | 0.3 | 3847 | 1691 | 5654 | 64 |
| 0.08 | 0.35 | 7904 | 995 | 2346 | 54 |
| 0.12 | 0.1 | 682 | 43164 | 55528 | 1074 |
| 0.12 | 0.15 | 715 | 30292 | 44026 | 1074 |
| 0.12 | 0.2 | 2203 | 14554 | 23210 | 386 |
| 0.12 | 0.25 | 5562 | 5871 | 10360 | 109 |
| 0.12 | 0.3 | 5916 | 3715 | 7662 | 101 |
| 0.12 | 0.35 | 12662 | 1944 | 2414 | 70 |

The error rate of the primer alignments in the first 1000 reads, with
0.2, as measured by `-e auto` (`auto-settings.c`): 0.9, 1.7, 3.0, 5.0 and
6.3% (below the mean of the simulation: the alignments with more errors
than allowed are not found, and the rate counts the read bases aligned).

The rule of `-e auto` (as in Lumipore: two values, one threshold): 0.2
below a measured rate of 4%, 0.3 from 4% on; it chooses 0.2 for the
reads at 1, 2 and 4%, and 0.3 at 8 and 12%. At 1-4%, 0.2 leaves 1073 to
1819 primer bases and 1 to 23 of the internal site, and removes 607 to
842 template bases; 0.3 would remove about 3 times more template for
about half the primer bases. At 8 and 12%, 0.2 leaves 6540 and 14554
primer bases and 157 and 386 of the internal site; 0.3 leaves 1691 and
3715, 64 and 101, and removes 3847 and 5916 template bases (1296 and
2203 with 0.2). 0.35 removes 2 to 3 times more template than 0.3 at
every rate. The consensus was not assembled for these runs.

## 2026-10-30: flucut with the adapter search

What was run: flucut at commit 6393a2d (adapters searched at the read
ends without a primer, `--adapters nb,rbk,lsk`, the default), `-j 1`,
on the reads of the sections below, scored with `score_reads.py`; on the
barcoded reads, IRMA 1.4.0 FLU-minion with and without the elongation,
compared with `compare_consensus.py`.

| reads | template removed | primer kept | flank kept | internal site kept |
|---|---|---|---|---|
| with barcodes, `--adapters none` | 254 | 2221 | 28539 | 99 |
| with barcodes, adapters | 267 | 1019 | 698 | 44 |
| without barcodes, `--adapters none` | 680 | 1819 | 7542 | 23 |
| without barcodes, adapters | 680 | 1819 | 7542 | 23 |

On the barcoded reads, an adapter was removed at the start of 234 reads
and at the end of 181. The reads without barcodes are trimmed as without
the search (no adapter found). Barbell, on the same barcoded reads:
1883, 445, 4375, 21; dorado trim: 3714, 6452, 360, 27 (below).

The IRMA consensus of the barcoded reads trimmed with the adapter search
is the same as before: with the elongation, 96 and 104 bases missing at
the ends (the primers), no error, ATA at the I348K codon; without it,
219 and 276 missing, no error.

Two settings were changed on the way, both found on these reads:
an adapter with a barcode was first allowed to overhang the read, and
then matched almost anywhere, the part left in the read being mostly
the barcode (Ns): on the reads without barcodes, 194 more template
bases were removed; and after the removal of an adapter, the primer was
first searched in the usual window of 100 bases, where primer-like
sequences of the template were taken for it (1255 template bases
removed): the primer must now lie within 5 bases of the adapter.

Time (3200 reads, one thread): 0.75 s on the barcoded reads, with or
without the search.

## 2026-10-30: flucut, one alignment for the primers at the read ends

What was run: flucut at commit 9e69ba0 (the search of a primer with
missing outer bases in one alignment instead of one per shorter form),
`-j 1`, on the reads of the next two sections, scored with
`score_reads.py`. IRMA was not run again.

| reads | template removed | primer kept | flank kept | internal site kept |
|---|---|---|---|---|
| without barcodes, before | 882 | 1852 | 7557 | 23 |
| without barcodes, 9e69ba0 | 680 | 1819 | 7542 | 23 |
| with barcodes, before | 251 | 2198 | 28685 | 110 |
| with barcodes, 9e69ba0 | 254 | 2221 | 28539 | 99 |

Time for the 3200 reads without barcodes, one thread: 0.63-0.71 s
(three runs, load average about 2.5), instead of 1.15-1.35 s.

## 2026-10-30: all the tools (simulated reads with barcodes)

What was run: `benchmark/run-tools.sh` as in the next section, on 3200
reads simulated with `benchmark/simulate.py --barcodes 4` (same genome
and seed defaults as the reads without barcodes; each amplicon carries
at both ends the native barcoding adapter of Oxford Nanopore of one of
NB01-NB04, outside the primers). The flank column counts the adapter,
barcode and random remnant bases. Barbell 0.3.3 (commit 057a902) was
run here, as it needs barcoded reads (see "Barbell" below); the other
versions are those of the next section.

Reads (bases by their true origin, `score_reads.py`):

| tool | written | lost | template removed | primer kept | flank kept | internal site kept |
|---|---|---|---|---|---|---|
| untrimmed | 3200 | 0 | 0 | 102899 | 423867 | 2731 |
| flucut | 3177 | 0 | 251 | 2198 | 28685 | 110 |
| cutadapt | 3184 | 0 | 107913 | 3371 | 29700 | 35 |
| IRMA-core, window 40 | 3200 | 0 | 0 | 102899 | 423785 | 2731 |
| IRMA-core, window 70 | 3199 | 0 | 23 | 80547 | 343053 | 2027 |
| IRMA-core, window 120 | 3200 | 0 | 121 | 13678 | 78267 | 705 |
| BBDuk, hdist 1 | 3200 | 0 | 97 | 17217 | 82495 | 717 |
| BBDuk, hdist 2 | 3200 | 0 | 360 | 12193 | 59862 | 557 |
| Porechop_ABI | 3186 | 7 | 18100 | 48424 | 212979 | 1298 |
| Porechop_ABI, tuned | 3179 | 7 | 30614 | 15361 | 73752 | 693 |
| dorado trim | 3200 | 0 | 3714 | 6452 | 360 | 27 |
| samtools ampliconclip | 3174 | 0 | 159 | 40 | 88 | 1337 |
| Barbell | 3174 | 1 | 1883 | 445 | 4375 | 21 |

(Porechop_ABI split 20 reads at adapters or primers found in their
middle.)

IRMA consensus with the reference elongation, summed over the 8
segments (columns as in the next section). Without the elongation every
tool gives, as without barcodes, 219 and 276 bases missing and no
error.

| tool | missing 5' | missing 3' | beyond the ends or errors | ambiguous | codon |
|---|---|---|---|---|---|
| untrimmed | 0 | 0 | 940 | 2 | ATA |
| flucut | 96 | 104 | 0 | 0 | ATA |
| cutadapt | 96 | 104 | 0 | 0 | ATA |
| IRMA-core, window 40 | 0 | 0 | 932 | 2 | ATA |
| IRMA-core, window 70 | 4 | 0 | 940 | 1 | ATA |
| IRMA-core, window 120 | 72 | 103 | 66 | 0 | ATA |
| BBDuk, hdist 1 | 72 | 58 | 95 | 0 | ATA |
| BBDuk, hdist 2 | 94 | 79 | 0 | 0 | ATA |
| Porechop_ABI | 14 | 60 | 646 | 1 | ATA |
| Porechop_ABI, tuned | 46 | 103 | 229 | 1 | ATA |
| dorado trim | 95 | 104 | 0 | 0 | ATA |
| samtools ampliconclip | 96 | 104 | 0 | 0 | ATA |
| Barbell | 96 | 104 | 0 | 0 | ATA |

What this shows, on these reads:

- **The consensus.** With the elongation, flucut, cutadapt, dorado
  trim, samtools ampliconclip and Barbell give a consensus without
  errors that ends where the primers began (dorado 1 base short of
  that). BBDuk with hdist 2 has no error either, but ends 2 and 25
  bases inside the primer positions. The tools which leave the adapters
  in many reads (IRMA-core, BBDuk with hdist 1, Porechop_ABI) give
  66 to 940 bases beyond the ends or errors.
- **The I348K site.** Here the untrimmed reads give the right codon
  (ATA), unlike the reads without barcodes, but only by a few reads: at
  the middle base of the codon, 158 reads have T and 152 A (without
  barcodes: 137 T, 224 A). The minimap2 alignment of IRMA soft-clips
  the start of the internally primed reads, 67 bases on average after
  the adapter instead of 20 without it, and with it a part of the
  primer bases. After flucut: 157 T, 1 A. The artefact thus depends on
  the library and on chance; the trimmed reads do not.
- **The reads.** flucut removes 251 template bases and leaves 2198
  primer bases; only the k-mer tools (IRMA-core, BBDuk with hdist 1),
  which leave 6 to 47 times more primer bases, and ampliconclip remove
  fewer template bases. Barbell, made for barcoded reads and given the adapter,
  barcode and primer as one query, leaves the fewest primer (445) and
  flank (4375) bases after ampliconclip, at the cost of 1883 template
  bases; dorado, with the kit, removes the adapters best (360 flank
  bases left). flucut leaves 28685 flank bases, at the read ends where it
  finds no primer (the flank before a primer is removed with it):
  flucut does not search the adapters.
  These reads normally reach flucut already demultiplexed and trimmed of
  the barcodes, e.g. by dorado or MinKNOW.

## 2026-10-30: all the tools (simulated reads, without barcodes)

What was run: `benchmark/run-tools.sh` (the commands of each tool are
in the script) on the 3200 simulated reads of the IRMA checks; flucut
at commit 3bfb35c, cutadapt 5.2, IRMA-core 0.11.0 (trimmer), BBTools
40.05 (BBDuk), Porechop_ABI 0.5.1 (commit 0bc9f17), dorado 2.1.2,
samtools 1.21 and minimap2 2.29 (those of IRMA), IRMA 1.4.0
(FLU-minion); the programs are given to the script in environment
variables.

Reads (bases by their true origin, `score_reads.py`):

| tool | written | lost | template removed | primer kept | flank kept | internal site kept |
|---|---|---|---|---|---|---|
| untrimmed | 3200 | 0 | 0 | 102554 | 127032 | 2698 |
| flucut | 3176 | 0 | 882 | 1852 | 7557 | 23 |
| cutadapt | 3177 | 0 | 102408 | 3037 | 8784 | 25 |
| IRMA-core, window 40 | 3195 | 0 | 77 | 43474 | 75068 | 1105 |
| IRMA-core, window 70 | 3200 | 0 | 130 | 13152 | 22890 | 660 |
| IRMA-core, window 120 | 3200 | 0 | 545 | 13133 | 22914 | 660 |
| BBDuk, hdist 1 | 3194 | 0 | 261 | 12974 | 22826 | 660 |
| BBDuk, hdist 2 | 3194 | 0 | 2977 | 7814 | 15831 | 502 |
| Porechop_ABI | 3182 | 6 | 28349 | 46964 | 62672 | 1389 |
| Porechop_ABI, tuned | 3179 | 6 | 81812 | 9107 | 14313 | 25 |
| dorado trim | 3200 | 0 | 2364 | 4904 | 11772 | 41 |
| samtools ampliconclip | 3173 | 0 | 142 | 30 | 135 | 1549 |

(Porechop_ABI split 11 reads at primers found in their middle.)

IRMA consensus, summed over the 8 segments, with the reference
elongation (`SKIP_E=0`); the codon of the I348K site in the majority
consensus. Without the elongation, every tool gives the same ends (219
and 276 bases missing) and no error. With the elongation, the bases of
the consensus beyond the 5' end of a segment are counted by the
alignment either as bases beyond the end or as errors (equal cost):
the two columns are summed here.

| tool | missing 5' | missing 3' | beyond the ends or errors | ambiguous | codon |
|---|---|---|---|---|---|
| untrimmed | 0 | 0 | 358 | 3 | **AAA** |
| flucut | 96 | 103 | 0 | 0 | ATA |
| cutadapt | 96 | 103 | 0 | 0 | ATA |
| IRMA-core, window 40 | 0 | 0 | 365 | 1 | ATA |
| IRMA-core, window 70 | 84 | 101 | 11 | 0 | ATA |
| IRMA-core, window 120 | 84 | 101 | 11 | 0 | ATA |
| BBDuk, hdist 1 | 84 | 101 | 19 | 0 | ATA |
| BBDuk, hdist 2 | 96 | 101 | 0 | 0 | ATA |
| Porechop_ABI | 14 | 59 | 228 | 1 | ATA |
| Porechop_ABI, tuned | 96 | 103 | 0 | 0 | ATA |
| dorado trim | 96 | 103 | 0 | 0 | ATA |
| samtools ampliconclip | 96 | 103 | 0 | 0 | ATA |

(96 and 103: the 12 and 13 bases of the primers in each segment; 103
and not 104, as an N appended by IRMA to segment 3 matches any base.)

What this shows, on these reads:

- **The consensus.** Every trimming tool removes the I348K artefact.
  With the elongation, six tools give the same consensus as flucut,
  without errors and ending where the primers began: cutadapt, BBDuk
  with hdist 2, Porechop_ABI tuned, dorado trim and samtools
  ampliconclip. The tools which leave primers in many reads (IRMA-core,
  BBDuk with hdist 1) leave 11 to 19 errors at the ends; IRMA-core with
  window 40 and Porechop_ABI with its defaults let the consensus extend
  beyond the segments.
- **The reads.** Among the tools which search the primer sequences,
  flucut leaves the fewest primer and flank bases, with fewer template
  bases removed than all but the k-mer tools (IRMA-core, BBDuk with
  hdist 1), which leave 7 to 23 times more primer bases. samtools
  ampliconclip is the best at the ends (30 primer bases left), with the
  reads mapped to the exact genome they were simulated from (a real
  sample differs from its reference); it leaves 1549 of the 2698 bases
  of the internal priming site, as its primer intervals are at the
  segment ends only.

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

The other tools of `TOOLS.md`, other shares of internal priming, other error rates,
and real data.
