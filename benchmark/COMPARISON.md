# flucut versus other trimmers

What flucut does that the other primer/adapter trimmers do not, and where
it does not. Every number here is from `benchmark/RESULTS.md` (the
commands, versions and data are there); nothing is asserted without a
measurement, and the ties and the places where flucut loses are stated.

## The short answer

For the **consensus** at the internal-priming rates seen in the public
datasets measured so far (about 2 internally primed reads per full-length
read at the HA site), every primer-aware trimmer removes the Uni12
footprint and IRMA then calls the correct codon: flucut, cutadapt, BBDuk,
the IRMA-core trimmer and dorado all turn the I348K artefact (AAA, Lys)
into the true ATA (Ile). On that outcome there is **no advantage**: if a
pipeline only needs the consensus and already runs one of these tools, it
gets the same answer.

flucut's advantages are at the **read level** (how cleanly it trims,
which matters at low coverage and for minor-variant frequencies), in
**what it reports** (it detects and quantifies internal priming; the
others trim silently), and in being **specific to the influenza protocol
with no per-run configuration**. The read-level cleanliness translates
into a wrong consensus for the other tools only under heavier internal
priming than these datasets show; that boundary is quantified below
(derived from the measured retention; an IRMA validation at controlled
ratios is pending).

## Read-level cleanliness (real Nanopore run SRR36760543)

Two things decide a trimmer's quality at the internal-priming site: how
many **genuine** reads (from full-length amplicons, which show the true
ATA) it keeps, and how many **primed** reads (which carry the primer's
AAA) it removes. On the HA codon of this run (869 genuine, 1720 primed):

| tool | genuine kept (of 869) | primed removed (of 1720) |
|---|---|---|
| **flucut** | **867** | **1713** (7 left) |
| cutadapt (`-g`) | 632 | 1715 (5 left) |
| BBDuk | 865 | 1594 (126 left) |
| IRMA-core | 869 | 1554 (166 left) |
| dorado trim | 867 | 1717 (3 left) |

flucut and dorado are the only tools clean on **both** axes. cutadapt
loses 27% of the genuine reads; BBDuk and IRMA-core leave 7-10% of the
primed footprints in. dorado matches flucut but is far slower and heavier
(below) and needs a primer file.

**Flip point.** From these retention rates, the consensus turns to the
wrong AAA when the primed:genuine ratio exceeds roughly: flucut 245,
cutadapt 250, BBDuk 14, IRMA-core 10. This run is at 2.0, so all agree.
flucut thus tolerates about 18x more internal priming than IRMA-core, and
13x more than BBDuk, before the consensus breaks. (Derived from the
measured retention; validation with IRMA at 2x/12x/20x is pending.)

## Per tool

### cutadapt

- The closest competitor: fast, light, removes the footprint cleanly.
- **flucut's advantage:** cutadapt cannot keep the genuine reads *and*
  remove the internal footprints at once. Unanchored (`-g`) it cuts the
  Uni12-like sequence inside full-length reads and loses real template
  (632/869 Nanopore, 25/35 Illumina genuine reads kept); anchored
  (`-g ^`) it keeps them but then leaves the internal footprints in
  (25 of ~223 left, Illumina). flucut does both (867/869 and 7 left; or
  35/35 and 3 left) because it searches the primer only in a window at
  the read end and uses the reference to tell a terminal start from an
  internal one. On simulated reads cutadapt removes 102408 template bases
  (its 5' primer matches internally) against flucut's 882.
- Also flucut-only: internal-priming report, complete segment ends with
  references (below), built-in scheme (no `-g/-a/-e/--revcomp/-n` to
  tune), one IRMA-ready output.
- **Where cutadapt is as good or better:** the consensus at ordinary
  priming rates (a tie), and speed (3.1 s vs 7.4 s on the Nanopore run;
  flucut's time includes the reference classification).

### BBDuk

- k-mer primer matching at the read ends, two passes (ktrim=l, ktrim=r).
- **flucut's advantage:** BBDuk under-trims the footprint (126/1720 left
  on the real Nanopore run, so its flip point is ~14x against flucut's
  ~245x), and it uses far more memory (about 950 MB against flucut's
  60 MB, from the JVM). On simulated reads BBDuk with hdist 2 leaves
  7814 primer bases against flucut's 1852. Plus the flucut-only report,
  complete ends with references, and the built-in scheme.
- **Where BBDuk is as good or better:** the consensus at ordinary rates
  (a tie), and wall-clock speed (3.9 s against flucut's zlib-era 4.0 s;
  flucut with libdeflate is now 1.6 s).

### IRMA-core trimmer

- The trimmer bundled with IRMA-core: k-mer primers, a restrict window.
- **flucut's advantage:** it under-trims the most of the lightweight
  tools (166/1720 left on the real Nanopore run, flip point ~10x, the
  lowest here), so it is the first to give a wrong consensus as internal
  priming rises. flucut also reports the priming and can keep the
  termini.
- **Where IRMA-core is as good or better:** the consensus at ordinary
  rates (a tie); it keeps all genuine reads (869/869); and it uses the
  least memory (11 MB) as it is single-threaded and streaming.

### dorado trim

- Oxford Nanopore's own trimmer, primers given as a primer file.
- **flucut's advantage:** speed and memory. On 32000 reads dorado trim
  takes 8 s and 480 MB against flucut's 2.4 s (zlib) / 1.6 s (libdeflate)
  and 60 MB; on the real Nanopore run 46 s against 7 s. flucut also
  reports the internal priming, keeps the termini on request, and needs
  no kit or primer file. On accuracy at the site the two are equal
  (dorado 867 genuine / 3 primed left).
- **Where dorado is as good or better:** read-level accuracy is equal,
  and dorado demultiplexes (flucut does not).

### Barbell

- A general Nanopore barcode/adapter finder; it also demultiplexes.
- **flucut's advantage:** on barcoded simulated reads Barbell removes
  1883 template bases against flucut's 251 (with adapter search, 267);
  it needs barcoded reads and a per-kit query and filter configuration
  (three steps: annotate, filter, trim), while flucut runs in one command
  on reads with or without barcodes; and flucut reports internal priming
  and can keep the termini. The consensus is the same for both.
- **Where Barbell is as good or better:** it leaves fewer primer (445 vs
  flucut 1019 with adapter search) and adapter bases, it demultiplexes
  (and on the simulated barcoded reads assigns 3174/3200 with no errors,
  against dorado's 3200/3200), and it is as fast (2.4 s). For getting
  barcodes off reads that are not yet demultiplexed, Barbell is the
  better tool; flucut expects reads already demultiplexed (by MinKNOW,
  dorado or Barbell).

### Porechop_ABI

- An adapter trimmer, primers given as custom adapters.
- **flucut's advantage:** much faster and lighter (Porechop_ABI tuned
  took 34 s and 250 MB on 32000 reads against flucut's 2.4 s and 60 MB),
  it does not split reads at spurious internal matches (Porechop split
  several), and it reports the priming and keeps the termini. Untuned,
  Porechop leaves most primers in; tuned, it removes a large amount of
  template (81812 bases on the simulated reads).
- **Where Porechop is as good or better:** nothing measured; it is a
  general adapter tool, not built for this.

### samtools ampliconclip (mapping-based)

- Maps the reads and clips primer intervals from a BED of the segment
  ends.
- **flucut's advantage:** ampliconclip clips only the fixed BED intervals
  at the segment termini, so it does **not** remove the internal-priming
  footprint at all (it left 1549 of 2698 internal-site bases on the
  simulated reads, 1337 on the barcoded ones): it cannot correct the
  artefact. It also needs the exact genome, a mapping step, and a BED.
- **Where ampliconclip is as good or better:** at the true termini it is
  the most precise (30 primer bases left), when the reads map to the
  exact reference.

## Features no other tool here has

- **Detection and reporting of internal priming:** with references,
  flucut classifies each read end as terminal or internally primed,
  reports the sites, their reads and the segments, draws the coverage of
  each segment by read origin, and warns about amplicons shorter than
  their reference. The others trim silently, so a user cannot see that
  the artefact was present or how strong it was. This needs
  **complete-segment references, with the UTRs** (not coding sequences);
  flucut warns when the references look like coding sequences. On three
  real Nanopore runs it measured 12-46% of the HA 5' ends as internally
  primed, at the I348K site (`RESULTS.md`).
- **Complete segment ends with references (the default when `-r` is
  given):** keep the conserved 12/13-base termini at the ends
  confidently classified as terminal, while removing the whole primer at
  internal or unclassified ends ("cut when unsure", so the artefact
  never survives). This gives complete ends and artefact removal in one
  pass. Only flucut knows, inside a primer, where the tail ends and the
  conserved terminus begins, so only flucut can make this cut; a
  sequence matcher removes all or nothing.
- **Built for the protocol:** the Uni12/Uni13 primers and schemes are
  built in; one command, no primer sequences or per-kit parameters.
- **One IRMA-ready output per sample**, an HTML report and a JSON
  summary, SAM-style tags on each read, and a live mode that follows a
  MinKNOW run while it sequences.

## Speed and memory (32000 simulated reads, 4 cores)

Measured together on one machine (`benchmark/RESULTS.md`). flucut here is
the zlib build; with libdeflate it is now 1.6 s / 4 threads.

| tool | wall s | peak MB |
|---|---|---|
| flucut, 4 threads | 2.4 (1.6 with libdeflate) | 60 |
| cutadapt | 3.7 | 41 |
| BBDuk | 4.0 | 954 |
| IRMA-core (1 thread) | 7.5 | 11 |
| dorado trim | 8.1 | 479 |
| minimap2 + ampliconclip | 15.0 | 125 |
| Porechop_ABI tuned | 33.6 | 253 |
| Barbell (barcoded reads) | 2.4 | 44 |

## Caveats

- Most numbers are from simulated reads; the real-data numbers so far are
  8 public Illumina runs and 3 Nanopore runs of H5N1, all Uni12/Uni13
  amplicons (`RESULTS.md`). The flip points above are derived from the
  measured read retention, not yet observed to change a real consensus:
  on the heaviest real run so far (SRR36760520, 46% of the HA 5' ends
  internally primed) every primer-aware tool still calls the correct
  codon, well below flucut's ~99.6% flip point.
- The read-level advantages matter for the consensus only at low coverage
  or high internal priming, and for minor-variant frequencies; at
  ordinary coverage and priming the consensus is a tie with the other
  primer-aware tools.
