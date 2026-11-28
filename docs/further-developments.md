# Possible further developments

Rewritten after each task. The order is roughly the order in which the
work is planned.

The focus is the trimming of the primers and the preparation of the
reads for IRMA: an advantage over the other trimmers has to show there,
in the IRMA consensus, and is measured by the benchmark. The reports
support the trimming; they are not the purpose of the program.

## Next

- **Influenza B on a sequenced sample**: done on four public Illumina
  runs of B/Victoria (BioProject PRJNA1426172, paired reads of
  fragmented amplicons of the IBV-GA2 primers;
  `../influenza_b_real_runs/` of the manuscript); still to do is a
  Nanopore run, which would show how often the primer is not found at
  an end of a whole amplicon.
- **flucut-posthoc and influenza B**: the post-hoc check (`posthoc/`)
  is still that of influenza A: it searches the two conserved termini
  of influenza A next to the boundary of each alignment, numbers the
  segments as influenza A does and reads the subtype from the
  descriptions. Influenza B needs its own design for the search: the
  bases its primers add after the universal terminus are the sequence
  of the segment, so they align to the reference as template does, and
  the boundary of the alignment is not where the terminus ends.

- **The letters of the test groups**: `fc-tests` gives each group a
  letter, and three of them were letters of the common options (`-q`
  quiet, `-V` version, `-L` log file) while `-c` named two groups. The
  letters are now free and checked at every run; `make test_units`
  (`fc-tests -q`) had been running one group of the thirty. Worth a
  look for what else the silence hid.

- **The kept terminus against the reference**: at an end classified as
  terminal, the bases kept are taken from the primer alignment alone.
  Where a sibling primer of a cocktail matches better than the one
  which primed (the right-end primers of the PB1/PB2 pair and of the
  PA pair of influenza B differ by one base), the bases kept can be
  those of the other primer. Comparing the kept bases with the
  reference, as the extension of the alignment described under "Reads
  cut too deep" would, would catch it.

- **Real runs, next**: three real H5N1 Nanopore runs and one Illumina run
  are done (`benchmark/RESULTS.md`); still to do are the other 7 Illumina
  SRA runs whose GenBank HA consensus has the I351K codon, with
  `run-real-illumina.sh` generalised.
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
  5' primer matches inside the reads, e.g. at the I351K site); the
  IRMA-core trimmer leaves many primers. To measure where the template
  matters: low coverage, variant frequencies.

- **`-e auto`, next**: calibrated on simulated reads only (two values,
  one threshold); to check on real runs of R9 and R10 flow cells, and
  whether the error rate of the reads should also set the window or the
  probe errors.
- **A template base taken for the last primer base**: where the read
  cannot tell primer from template at the inner end of a primer (ends of
  equal cost), the base is kept as template. One template base is still
  removed where the read lost a primer base equal to the template base
  after it (one G deleted where the primer's last G meets a template G,
  as at the I351K site): the alignment with fewer errors then takes the
  template G for the primer's. On the simulated reads of the benchmark,
  12 bases in 3200 reads at 4% errors, 79 at 12%
  (`benchmark/RESULTS.md`, 2026-11-28). With references, the template
  base after the primer is known, and such a shortened run could be
  given to the primer. Not done: few bases.
- **Reads cut too deep** (6 of 3200 reads, about 70 template bases each):
  internally primed reads whose primer was cut to 4-5 bases after a
  flank; a primer-like sequence of the template within the window is
  taken for the primer. With references, an extension of the alignment
  to the reference from the read end would recognise it, and would also
  cut the ends without a primer found at the segment terminus (1864
  primer bases left in about 210 of 6400 read ends; the consensus is not
  affected in the benchmark). Not done: few reads.
- **Dashboard and report**: the dashboard (`flucut_dashboard.html`)
  follows the run and stays at its end; the report is written at the
  end. Both show the run overview (a row per sample, a cell per
  segment), the HA in all samples, and each sample (state, subtype,
  reads by origin, segments, read lengths and qualities, sequencing
  metadata). Still possible: a test on a replay of a real run; the yield over time from the start times of the
  reads; the remnants (primer bases
  left) and the primers found whole or partial; the thresholds of the
  states (depth 50, 90% of the segment) as options.
- **Reads of a second HA or NA subtype, next**: the report gives them
  (`src/analysis/subtype-reads.h`), counted in the reads the references are
  chosen from: the first 2000 reads of the first file of a sample,
  which in the files MinKNOW writes are often fewer. Still possible:
  counting them over all the reads of a sample, which would see smaller
  fractions and would not depend on its first file; a view of the run
  which shows which samples carry each other's subtypes, the pattern
  barcode crosstalk leaves; a measure on sequenced mixtures of known
  composition (the threshold comes from simulated reads and from the
  background of two sequenced runs); short reads, whose probes overlap
  and so vote alike, counted with fewer probes; a divergent NA, whose
  probes scatter over several subtypes (two H5N8 samples gave a remark
  for the NA while only a fifth of their NA reads agreed on N8). When
  the reference of the HA or of the NA is chosen, the votes are still
  counted together by the whole subtype name of the references (an
  H5N1 and an H5N2 reference of the HA apart), not by H or N number as
  here: to be made the same, with the subtype naming measured again.
- **Read-origin plot per segment**: done. For each reference, a stacked
  coverage profile along its length, built from flucut's classification
  (the span of each located read, from its probe position and kept
  length), coloured by origin (terminal, internal, by strand) with the
  reads of unknown 5' origin at the bottom; internal priming shows as a
  step up at the site, marked. It needs complete-segment references
  (with the UTRs); a warning fires when they look like coding sequences.
  Still possible: a mark of the codon of the site on the plot, and a
  separate histogram of the amplicon start positions.

## Benchmark

- Simulated reads with a known truth (all eight segments, several
  subtypes, the primer variants, tails and adapters, partial primers,
  dimers, chimeras, reads without primer; Nanopore error profiles).
- The tools of `benchmark/TOOLS.md`, each in its best documented
  configuration; **Barbell** is the reference to match or beat, and its
  `--use-extended` mode has to be included.
- Real data: public runs amplified with the universal primers (done: three
  H5N1 Nanopore runs and one Illumina run, `benchmark/RESULTS.md`).
- The measures which matter for adoption: primer and adapter bases left
  in the reads, template bases wrongly removed, reads lost, and the
  effect on the IRMA consensus at the segment ends. Speed and memory are
  secondary.

## Trimming and detection

- **Is it only HA?** (influenza A) Measured on real data: three H5N1 Nanopore runs and
  one Illumina run (`benchmark/RESULTS.md`) were classified against the
  references of all eight segments. The internal priming is concentrated
  in the HA, at the Uni12-like site (nt 1070-1094, the I351K codon), the
  same site in all three Nanopore runs; the other segments show only a few
  percent. So on this material it is essentially HA, and the simulation
  (`benchmark/simulate.py`, HA-only internal site) matches. Still possible:
  a per-segment, per-primer (Uni12 and Uni13) table to back the claim with
  numbers rather than prose, and other subtypes (so far only H5N1).
- **Choice of the references for a run**: done. A built-in set of
  complete segments is the default, downloaded on first use and cached;
  each barcode is classified against the reference of each segment
  closest to its reads (chosen from a sample of them, segment by
  segment). It covers the HA of every subtype (H1 to H16), the NA of
  every subtype (N1 to N9) and the other six segments of six subtypes,
  and every reference begins and ends with the conserved termini
  (`scripts/find_builtin_references.py` searches NCBI for such records).
  Still possible: a packaged copy for offline installs, a check that the
  unclassified ends stay rare, and estimating the position of the
  terminus in each reference from the reads, so that a reference without
  its UTRs (the seed sequences of IRMA, say) can be used as well.
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
  binding sites on the segments. A user-defined primer would also need
  the segments it binds, as the primers of influenza B have.

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
