# Possible further developments

Rewritten after each task. The order is roughly the order in which the
work is planned.

The focus is the trimming of the primers and the preparation of the
reads for IRMA: an advantage over the other trimmers has to show there,
in the IRMA consensus, and is measured by the benchmark. The reports
support the trimming; they are not the purpose of the program.

## Next

- **Real runs, next**: the other 7 SRA runs whose GenBank HA consensus
  has the I348K codon (benchmark/RESULTS.md), with
  `run-real-illumina.sh` generalised; a real Nanopore run.
- **Paired Illumina reads**: R1 and R2 are trimmed as two inputs; a
  read too short is dropped from one file only (IRMA accepted it on
  SRR34913296); keep the pairs in sync (an empty read instead).
- **Benchmark, next** (`benchmark/RESULTS.md`): other error rates and
  shares of internal priming, several simulations per setting (with
  barcodes, the artefact of the untrimmed reads came down to 158
  reads against 152); the frequencies of the variants; real data.
- **Where flucut trims better than a generic trimmer** (first results
  in `benchmark/RESULTS.md`): the consensus of cutadapt is the same on
  the simulated reads; cutadapt removes 80 times as much template (its
  5' primer matches inside the reads, e.g. at the I348K site); the
  IRMA-core trimmer leaves many primers. To measure where the template
  matters: low coverage, variant frequencies.

- **`-e auto`, next**: calibrated on simulated reads only (two values,
  one threshold); to check on real runs of R9 and R10 flow cells, and
  whether the error rate of the reads should also set the window or the
  probe errors.
- **Errors next to a primer** (measured, kept as it is): at equal cost,
  the alignment of a primer is extended inward, and can take one
  template base (e.g. the last primer base deleted in the read). Without
  the extension, on the reads of the benchmark, 208 fewer template bases
  were removed but 317 more primer bases kept (also at 1% and 12%
  errors: `benchmark/RESULTS.md`): no gain.
- **Reads cut too deep** (6 of 3200 reads, about 70 template bases each):
  internally primed reads whose primer was cut to 4-5 bases after a
  flank; a primer-like sequence of the template within the window is
  taken for the primer. With references, an extension of the alignment
  to the reference from the read end would recognise it, and would also
  cut the ends without a primer found at the segment terminus (1864
  primer bases left in about 210 of 6400 read ends; the consensus is not
  affected in the benchmark). Not done: few reads.
- **Live mode, next**: the dashboard of Lumipore (a report which fits
  the screen without scrolling while the run is in progress; today the
  report itself reloads every 30 seconds); a test on a replay of a real
  run (Lumipore's `tests/fake_minknow`).
  The report model (`report.h`: refresh, running, dashboard, the
  settings decided automatically) and the worker threads are ready for
  it; the dashboard itself is still to be written.
- **Report, next**: the remnants (primer bases left), the read lengths
  before and after the trimming, the primers found whole or partial,
  per segment; a mark of the internal sites in the length histograms.

## Benchmark

- Simulated reads with a known truth (all eight segments, several
  subtypes, the primer variants, tails and adapters, partial primers,
  dimers, chimeras, reads without primer; Nanopore error profiles).
- The tools of `benchmark/TOOLS.md`, each in its best documented
  configuration; **Barbell** is the reference to match or beat, and its
  `--use-extended` mode has to be included.
- Real data: one or two public Nanopore runs amplified with the
  universal primers.
- The measures which matter for adoption: primer and adapter bases left
  in the reads, template bases wrongly removed, reads lost, and the
  effect on the IRMA consensus at the segment ends. Speed and memory are
  secondary.

## Trimming and detection

- **Is it only HA?** Run the classification on real data against the
  references of all eight segments and report the internal sites found
  per segment and per primer (Uni12 and Uni13); extend the defaults if
  sites appear outside HA. Not yet measured.
- **Choice of the references for a run**: a default set from
  `scripts/fetch_references.py` covering the subtypes seen at the
  institute, and a check that the unclassified ends stay rare with it.
- **Start contexts**: group the reads by the bases which follow the 5'
  primer and rank the groups; a second reference-free cross-check, next
  to the length check (`lengths.tsv`), for the report.
- **Length check without references**: the expected length of a
  segment (today only the peaks are listed, in the log; the report has
  one histogram of all lengths).
- **Concatemers**: several amplicons in one read, split into separate
  reads; chimeras (amplicons of different segments joined). The internal
  matches already found are their junctions.
- **Adapters, next**: the adapters are searched only at the ends
  without a primer, with the barcode as Ns (none identified); to do:
  real reads (the simulated adapters have no remnants between the
  barcode and the primer), the other kits of Dorado (PCR barcoding,
  16S), the barcode identified (demultiplexing, as Dorado and Barbell),
  and a check that reads with adapters in their middle (chimeras of the
  library) are counted.
- **Barcoded primers**: a barcode between the adapter and the primer,
  used to assign the reads to samples. The layouts of the published
  barcoded protocols still have to be read from their supplementary
  tables (see below).
- **User-defined schemes** (YAML): the tail, the adapter and the
  barcodes can be changed; the conserved termini cannot, as they are the
  binding sites on the segments.

## Illumina (second phase)

- Paired-end input, with the pairs kept in sync; the primers are only at
  the ends of the fragments which cover a segment end, as the fragments
  come from tagmentation.

## Documentation

- Manual (`docs/manual`, `manual/flucut-manual.pdf`): the results of
  all the tools of the benchmark in the performance chapter.
- Developer guide (`docs/dev_manual.md`).

## Open questions

- The exact ordered form of the degenerate Uni12 primer (the published
  degenerate form is built in; the supplementary table of the
  publication has still to be checked).
- The layouts of the published barcoded primer sets (supplementary
  tables).
- Whether the reads should optionally be written in a single
  orientation; at the moment they keep the orientation they were
  sequenced in, as a trimmer should not change it.

## Checks still to set up

- Fuzzing harnesses (FASTA reader, options), continuous integration
  started manually.
- The subtype filter of `scripts/fetch_references.py`: asked for H5N1
  with `--refseq-only`, it returned an H9N2 record (see
  `testdata/README.md`).
- A review of the code style and of the bugs, with subagents.
