# flucut versus other trimmers

What flucut does that the other primer and adapter trimmers do not, and
where it does not. Every number is from the entry of 2026-11-27 of
`benchmark/RESULTS.md` (flucut commit d26b9b1; versions and commands
are there) and from the tables it names, in
`/home/baia/prj/flucut/trimmer_benchmark_runs/tables/`. The ties and
the places where flucut loses are stated. dorado trim could not be run
for that entry: where it appears, the number is that of 2026-10-30 and
is marked so.

flucut appears in three configurations: with its default options (the
reads classified against the built-in references, the conserved termini
kept at the terminal ends), with `-r off` (no references, every primer
base removed, the primers of influenza A and B both searched) and with
`-r off --virus a` (the primers of influenza A alone).

The artefact is called I351K: an isoleucine read as lysine (ATA as
AAA) at residue 6 of HA2, the start of the fusion peptide, where the
Uni12 primer anneals inside the H5 HA; 351 is counted on the HA0
precursor from the initiator methionine.

## The short answer

For the **consensus**, every primer-aware trimmer removes the Uni12
footprint, and IRMA then calls the template codon: on the simulated
reads, on three public Nanopore runs in which 12 to 45% of the located
5' read ends are internally primed, and on a public Illumina run. On
that outcome there is **no advantage**: a pipeline which only needs
the codon and already runs one of these tools gets the same answer.

flucut's advantages are at the **read level** (how cleanly it trims),
in **what it reports** (it detects and measures internal priming; the
others trim silently), in the **complete segment ends** it gives with
references, and in being **specific to the protocol, with no
configuration**. At which share of internal priming the reads left by
the other tools would change a consensus has not been measured.

## Read-level cleanliness on public runs

At the internal priming site a trimmer should keep the reads of
amplicons started upstream, which show the template codon ATA, and
leave no read which starts at the site with the primer codon AAA.
Reads over the codon in three Nanopore runs
(`tables/public_site_reads.tsv`; template reads kept / reads with the
primer codon left):

| tool | SRR36760543 (869 / 1745) | SRR36760506 (1341 / 3009) | SRR36760520 (611 / 9208) |
|---|---|---|---|
| flucut, default | 870 / 0 | 1356 / 0 | 611 / 0 |
| flucut `-r off` | 869 / 36 | 1343 / 101 | 610 / 274 |
| cutadapt, not anchored | 632 / 6 | 932 / 21 | 482 / 30 |
| BBDuk | 865 / 147 | 1334 / 243 | 610 / 787 |
| IRMA-core trimmer | 869 / 190 | 1341 / 304 | 611 / 1033 |

(The reads with the primer codon are those starting at the site and
those starting elsewhere, summed. flucut can show a few more template
reads than the untrimmed reads, as a read cut at the right place maps
from further upstream.)

flucut with references is the only tool clean on both counts.
cutadapt, not anchored, loses 21 to 31% of the template reads, which
it cuts at the site. BBDuk and the IRMA-core trimmer leave 8 to 11% of
the reads with the primer codon. flucut without references leaves 2 to
3%.

## Per tool

### cutadapt

- The closest competitor: light, and it removes the footprint.
- **flucut's advantage:** cutadapt cannot keep the template reads and
  remove the internal footprints at once. Not anchored it cuts the
  Uni12-like sequence inside full-length reads (632 of 869, 932 of
  1341 and 482 of 611 template reads kept on the Nanopore runs, 25 of
  35 on the Illumina run; 8 to 11% of all bases of the Nanopore runs
  removed, against 2 to 4% by flucut); anchored, on the Illumina run,
  it keeps them and leaves 20 reads with the primer codon. On the
  simulated reads it removes 102408 template bases against flucut's
  680 to 1175.
- Also flucut only: the report of internal priming, complete segment
  ends, a built-in scheme, one output ready for IRMA.
- **Where cutadapt is as good or better:** the consensus (a tie); the
  primer bases left on the simulated reads when flucut runs with `-r
  off` (3037 against 2966, about the same; 1760 with `--virus a`).

### BBDuk

- k-mer matching of the primers at the read ends, two passes.
- **flucut's advantage:** BBDuk leaves the footprint in 8 to 9% of the
  primed reads of the Nanopore runs, and uses 400 to 470 MB of memory
  against 54 to 101 MB. On the simulated reads, with hdist 2, it leaves
  7814 primer bases and removes 2977 template bases.
- **Where BBDuk is as good or better:** the consensus (a tie), and the
  time on the simulated reads: 2.8 s against 4.4 s for flucut with its
  default options (1.3 s for flucut with `-r off --virus a`).

### IRMA-core trimmer

- The trimmer of IRMA-core: the primers as k-mers, within a window.
- **flucut's advantage:** it leaves the footprint in 10 to 11% of the
  primed reads of the Nanopore runs, the most of the tools compared,
  and 72 of 221 on the Illumina run. On the simulated reads it leaves
  13133 primer bases or more, and with the elongation of IRMA its
  consensus has errors at the ends.
- **Where the IRMA-core trimmer is as good or better:** the codon (a
  tie); it keeps all template reads; it uses the least memory (12 MB).

### dorado trim

- Not measured for the entry of 2026-11-27 (not installed). On
  2026-10-30, with dorado 2.1.2 on the simulated reads: a consensus
  without errors, 2364 template bases removed and 4904 primer bases
  left, 8.1 s and 479 MB for 32,000 reads. dorado also demultiplexes;
  flucut does not.

### Barbell

- A barcode and adapter finder for Nanopore reads; it demultiplexes.
- **flucut's advantage:** on the simulated barcoded reads Barbell
  removes 1883 template bases against 268 (flucut `-r off --virus a`);
  it needs barcoded reads and a query and filter configuration of
  three steps; flucut reports internal priming and can keep the
  termini. The consensus is the same.
- **Where Barbell is as good or better:** it leaves fewer primer bases
  (445 against 956) and it demultiplexes. For reads which are not yet
  demultiplexed Barbell is the tool; flucut expects demultiplexed
  reads.

### Porechop_ABI

- An adapter trimmer, the primers given as custom adapters.
- **flucut's advantage:** time and memory (25.7 s and 254 MB for
  32,000 reads); with its defaults Porechop_ABI leaves most primers
  in, and tuned it removes 81812 template bases of the simulated reads.
- **Where Porechop_ABI is as good or better:** nothing measured.

### samtools ampliconclip

- Maps the reads and clips primer intervals given in a BED file.
- **flucut's advantage:** the intervals are at the segment ends only,
  so the footprint of internal priming stays (1549 of 2698 bases of the
  internal site left on the simulated reads, 1337 of 2731 on the
  barcoded ones). It needs a reference, a mapping and a BED file.
- **Where ampliconclip is as good or better:** at the true ends it is
  the most precise (30 primer bases left, 142 template bases removed),
  given the exact genome the reads were simulated from. Its consensus
  has the template codon here, as the bases left are a minority.

## Features no other tool here has

- **Detection and report of internal priming.** With references (the
  default), flucut classifies each read end as terminal or internally
  primed, reports the sites and their reads, and draws the coverage of
  each segment by read origin. On three public Nanopore runs it found
  12.2, 13.1 and 45.1% of the located 5' ends internally primed, the
  largest sites in the HA around the I351K codon (63 to 92% of the HA
  reads are internally primed), and further sites in PB1, PB2 and PA. This needs
  complete segments as references, with the untranslated regions:
  against coding sequences the terminal ends are not located (91 of
  47,990 reads on SRR36760543), and flucut warns.
- **Complete segment ends.** The conserved terminus is kept at the
  ends classified as terminal and the whole primer removed elsewhere.
  With the elongation of IRMA the consensus of the simulated reads is
  then complete at both ends of all segments (every other tool lacks
  12 and 13 bases per segment or has errors there); its terminal bases
  are those of the primers (1 ambiguous base per segment, 2
  differences from the truth in 8 segments).
- **Reads with several amplicons are split** (173 simulated reads of
  two amplicons: 37 left whole without references, 4 with the genome as
  references).
- **Built for the protocol:** the primers and schemes are built in;
  one command, no primer file.
- **One output per sample ready for IRMA**, an HTML report and a JSON
  summary, tags on each read, and a live mode which follows a run.

## Where flucut's default configuration loses to its own `-r off --virus a`

On the simulated reads, whose genome (PR8 with an H5 HA) is far from
the built-in references in several segments:

- 1175 template bases removed instead of 680, and 60 of 3200 reads of
  one amplicon written in two pieces (0 without references, 0 with the
  genome as references); 153 and 183 at 8 and 12% errors;
- 4.4 s instead of 1.3 s for 32,000 reads, and 101 MB instead of 60.

On the public runs, whose viruses are close to a built-in reference,
the default is as fast as with 8 references. And `-r off` without
`--virus a` is the poorest of the three: the primers of influenza B
are searched in every read (4.5 times the time) and are sometimes
taken for the primer of an end (2966 primer bases left instead of
1760, 83 bases of the internal site instead of 23).

## Time and memory

32,000 simulated reads, each tool confined to 4 processors of a shared
64-core machine, median of three runs
(`tables/speed_32000_reads.tsv`):

| tool | wall s | peak MB |
|---|---|---|
| flucut `-r off --virus a`, 4 threads | 1.28 | 60 |
| flucut `-r` 8 references, 4 threads | 1.75 | 88 |
| BBDuk | 2.79 | 439 |
| cutadapt | 4.10 | 48 |
| flucut, default, 4 threads | 4.41 | 101 |
| flucut `-r off`, 4 threads | 5.77 | 63 |
| IRMA-core trimmer (1 thread) | 6.32 | 12 |
| minimap2 + ampliconclip | 12.13 | 124 |
| Porechop_ABI tuned | 25.72 | 254 |
| Barbell (barcoded reads) | 3.53 | 68 |

On the public Nanopore runs (48,000 to 77,000 reads;
`tables/public_speed.tsv`): flucut with its default options 2.2 to
3.4 s, BBDuk 2.9 to 3.3 s, cutadapt 3.5 to 5.1 s, the IRMA-core
trimmer 7.3 to 9.3 s. At the scale of a run (1,129,506 Nanopore reads,
4 threads): 38 s and 55 MB with the default options.

## Caveats

- Most numbers are from simulated reads of one constructed genome. The
  public data are three Nanopore runs of one project and one Illumina
  run, all H5N1.
- On every data set measured the consensus codon is the same for all
  primer-aware tools. The read-level differences matter at low
  coverage, at higher shares of internal priming and for the
  frequencies of minor variants; none of these was measured with the
  other tools beyond the low-coverage series, in which they tie.
- The times were measured on a shared machine, with other jobs running
  on other cores.
