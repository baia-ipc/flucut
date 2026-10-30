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
