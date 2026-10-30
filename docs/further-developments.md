# Possible further developments

Rewritten after each task. The order is roughly the order in which the
work is planned.

The focus is the trimming of the primers and the preparation of the
reads for IRMA: an advantage over the other trimmers has to show there,
in the IRMA consensus, and is measured by the benchmark. The reports
support the trimming; they are not the purpose of the program.

## Next

- **Benchmark, other tools** (`benchmark/`): the IRMA checks are done
  (`benchmark/RESULTS.md`: flucut removes the I348K artefact, IRMA
  accepts the tags, the consensus ends 12 and 13 bases short with the
  default trimming). Next: the read-level measures from
  `truth.tsv.gz`, then cutadapt, Barbell (default and
  `--use-extended`), the IRMA-core trimmer, Porechop_ABI, BBDuk, dorado
  trim and the alignment-based trimmers, on the same reads; other
  shares of internal priming and error rates.
- **Where flucut can trim better than a generic trimmer**, to be
  verified with the benchmark, not assumed: the 12-base degenerate
  Uni12 primer found partial and with Nanopore errors at the read end
  (remnants left by trimmers which need longer or exact adapters); no
  template removed before a primer-like site inside a read (e.g. a
  non-anchored 5' adapter in cutadapt removes everything before a match
  anywhere in the read); both orientations; no configuration.

- **Per-barcode outputs of a MinKNOW run directory**: today a directory
  is one sample; a run directory (`fastq_pass/barcodeNN/`) should give
  one output per barcode, without listing the barcodes.
- **`-e auto`**: the primer error rate estimated from the first reads.
- **Errors next to a primer**: with a substitution right next to a
  primer, the alignment can extend over it and remove one template base
  (seen in the program tests, read `errors_plus`). Small; to be weighed
  against the remnants left by a tighter alignment, in the benchmark.
- **Internal priming sites**: the counts of `summary.tsv` per end and
  origin exist; the sites themselves (reference and position, with
  their number of reads) still have to be collected for the report.
- **HTML report, always written** (self-contained, light and dark, no
  internet): run parameters, reads at each step, reads per sample and
  per segment, primers found per segment and end (whole, partial,
  missing), the internal priming sites with their share of the reads,
  remnants, read lengths before and after trimming.
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
  primer and rank the groups; a reference-free cross-check of the
  classification, for the report.
- **Concatemers**: several amplicons in one read, split into separate
  reads; chimeras (amplicons of different segments joined). The internal
  matches already found are their junctions.
- **Sequencing adapters** left by the basecaller, searched outside the
  primers.
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

- Sphinx manual, with the primers and their sources, the
  classification against the references, the trimming policies and
  what the consensus termini then contain.
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
  started manually, and the live mode tests (`make test_live`, not in
  `make test` until the live mode exists).
- The subtype filter of `scripts/fetch_references.py`: asked for H5N1
  with `--refseq-only`, it returned an H9N2 record (see
  `testdata/README.md`).
- A review of the code style and of the bugs, with subagents.
