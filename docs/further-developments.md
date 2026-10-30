# Possible further developments

Rewritten after each task. The order is roughly the order in which the
work is planned.

## Next

- **Start contexts** (the main measure of internal priming, see
  `docs/internal-priming.md`): group the reads of a run by the bases
  which follow the trimmed 5' primer and rank the groups, so that the
  starts of the real amplicons separate from the internal ones without
  any reference sequence; report the unexpected groups with their
  sequence and their share of the reads. Needs the counting pipeline.
- **Command line program** (`src/options.c`, `src/flucut.c`): inputs
  (FASTQ/BAM, files, directories, MinKNOW run directories), `--trim
  primer|terminus`, `-e auto` (errors from the first reads), minimum
  length after trimming, threads, output directory, the banner and the
  logo.
- **Output of the trimmed reads**: FASTQ (gzip) ready for IRMA, with
  what was found and removed at each end written in the read title as
  SAM-like tags; per sample/barcode when the input is a MinKNOW run
  directory.
- **Counting and statistics**: reads at each step, per end and per
  variant of the primer, primer and tail remnants, artefacts.
- **HTML report, always written** (self-contained, light and dark, no
  internet): run parameters, reads at each step, reads per sample and
  per segment, primers found per segment and end (whole, partial,
  missing), remnants, artefacts, read lengths before and after
  trimming.
- **JSON summary**, same content, versioned format, written to a
  temporary file and renamed.
- **Live mode** `--live`: follow a MinKNOW run while it is sequencing
  (`--min-age`, `--scan-interval`, `--resume`, end of the run from the
  final summary file, barcode directories appearing during the run);
  after each file the JSON summary is rewritten and the HTML report
  becomes a dashboard which fits the screen without scrolling and
  reloads itself every 1 to 30 seconds, and which becomes the final
  report at the end of the run. Test: a fast replay
  (`--scan-interval 0.05`) whose outputs must be identical to those of
  a normal run on the same reads.

## Trimming and detection

- **Concatemers**: several amplicons in one read, split into separate
  reads; chimeras (amplicons of different segments joined). The internal
  matches already found are their junctions.
- **Classification of the read starts against references** (optional,
  the most direct measure of internal priming when the segment
  references of the run are at hand).
- **Sequencing adapters** left by the basecaller, searched outside the
  primers.
- **Barcoded primers**: a barcode between the adapter and the primer,
  used to assign the reads to samples. The layouts of the published
  barcoded protocols still have to be read from their supplementary
  tables (see below).
- **User-defined schemes** (YAML): the tail, the adapter and the
  barcodes can be changed; the conserved termini cannot, as they are the
  binding sites on the segments.
- **Segment of a read**: which of the eight segments and which end a
  read covers, for the report (not needed for the trimming itself).
- **Padding of the termini**: with the default trimming the consensus is
  shorter than the segment at both ends; an option could write the
  termini as N, for the programs which expect complete ends (e.g. the
  FLU-utr module of IRMA). To be decided with the users.

## Illumina (second phase)

- Paired-end input, with the pairs kept in sync; the primers are only at
  the ends of the fragments which cover a segment end, as the fragments
  come from tagmentation.

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

## Documentation

- Sphinx manual, with the primers and their sources, the trimming
  policies and what the consensus termini then contain.
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

- Program tests (`testdata/cases.tsv`), option validation tests, fuzzing
  harnesses, continuous integration started manually.
- A review of the code style and of the bugs, with subagents.
