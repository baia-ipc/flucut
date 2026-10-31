# Figures of the user manual

- `make_diagrams.py`: the diagrams `read_layout` and `internal_priming`
  (PDF for the manual, PNG for the HTML version); requires matplotlib:

      python3 docs/manual/figures/make_diagrams.py

  The sequences are those of the primers (`src/primers.h`) and of the HA
  of KF918483.1 at the internal priming site (see
  `benchmark/README.md`).

- `title.jpg`: the illustration of the title page, provided by the
  author (an AI-generated image), converted from PNG to JPEG (quality
  90) with Pillow.

- `terminal.js`: a terminal transcript as an image (ANSI colors), from
  the manual of Lumipore, with the long lines wrapped; requires
  Playwright with Chromium and the fonts JetBrains Mono and FreeMono
  (Braille). `terminal_flucut.png`, from a copy of `testdata` in an
  empty directory, with `bin` in the PATH:

      script -qec "stty cols 100; flucut -o out -r testdata/ref-ha.fasta \
        testdata/reads.fastq" /dev/null > flucut.ans
      node terminal.js flucut.ans terminal_flucut.png \
        "flucut -o out -r testdata/ref-ha.fasta testdata/reads.fastq" 100

- `screenshot.js`: screenshots of the sections of an HTML report
  (`report_*.png`), from the manual of Lumipore (unchanged, but for its
  header); requires Playwright with Chromium. The reports, from a
  directory with copies of the reads and the genome of the simulation
  and of `testdata`, with `bin` in the PATH:

      flucut -q -o out -r genome.fasta reads.fastq.gz
      flucut -q -o out2 testdata/fastq_pass
      node screenshot.js out/flucut_report.html . report_head=head \
        report_tiles="+All reads" report_steps="Reads at each step" \
        report_internal="Internal priming" \
        report_complete="+Complete amplicons" \
        report_lengths="Lengths of the complete amplicons"
      node screenshot.js out2/flucut_report.html . report_groups=Groups

  The reads (`reads.fastq.gz`, 3200 reads, 239 of them started at the
  internal site of the HA) and the genome (`genome.fasta`) are those of
  the IRMA checks of the benchmark (`benchmark/README.md`): the genome
  downloaded by `benchmark/make-genome.sh genome`, the reads simulated
  with the defaults of `benchmark/simulate.py` (seed 1, 400 reads per
  segment):

      benchmark/make-genome.sh genome
      python3 -I benchmark/simulate.py genome/genome.fasta sim

- `report_coverage.png`: the coverage-by-read-origin section, from a real
  Nanopore run (it is more telling than the simulated reads, whose
  amplicons are nearly all full length). The run is SRA `SRR36760543`
  (an H5N1 isolate); the references are complete H5N1 segments **with the
  UTRs** (`coverage_refs.fasta`, the eight segments of one isolate,
  accessions `QB222556`-`QB222563`, fetched from NCBI), as the coverage
  needs the 5' UTR to place the terminal reads (CDS-only references
  cannot, see the manual). From a directory with the run's reads:

      flucut -q -o cov -r coverage_refs.fasta SRR36760543.fastq.gz
      node screenshot.js cov/flucut_report.html . \
        report_coverage="+Coverage by read origin"

- `dashboard.png`: the live-mode dashboard, from a live run of two
  barcode directories (the reads of two real Nanopore runs, `SRR36760543`
  and `SRR36760506`, split into `barcode01` and `barcode02`), against the
  complete-segment references (`coverage_refs.fasta`). Start a live run,
  drop a file into each barcode directory, and screenshot at a wide
  viewport while it is running:

      flucut --live --min-age 0 -r coverage_refs.fasta -o out in &
      # ... in/barcode01/*.fastq.gz, in/barcode02/*.fastq.gz appear ...
      node dashboard-shot.js out/flucut_report.html dashboard.png

  (`dashboard-shot.js`: Playwright, a 1680x945 viewport, a full-page
  screenshot.)
