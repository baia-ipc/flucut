# Figures of the user manual

- `make_diagrams.py`: the diagrams `read_layout` and `internal_priming`
  (PDF for the manual, PNG for the HTML version); requires matplotlib:

      python3 docs/manual/figures/make_diagrams.py

  The sequences are those of the primers (`src/seq/primers.h`) and of the HA
  of KF918483.1 at the internal priming site (see
  `benchmark/README.md`).

- `title.jpg`: the illustration of the title page, provided by the
  author, converted from PNG to JPEG (quality 90) with Pillow.

- `terminal.js`: a terminal transcript as an image (ANSI colors), with
  the long lines wrapped; requires
  Playwright with Chromium and the fonts JetBrains Mono and FreeMono
  (Braille). `terminal_flucut.png`, from a copy of `testdata` in an
  empty directory, with `bin` in the PATH:

      script -qec "stty cols 100; flucut -o out -r testdata/ref-ha.fasta \
        testdata/reads.fastq" /dev/null > flucut.ans
      node terminal.js flucut.ans terminal_flucut.png \
        "flucut -o out -r testdata/ref-ha.fasta testdata/reads.fastq" 100

- `report_figures.py`: the figures of the HTML report and of the
  dashboard (`report_overview`, `report_ha`, `report_sample`,
  `report_segments`, `report_reads`, `dashboard`): it renders the pages
  with a headless Chromium and cuts them at their headings; requires
  Chromium (the program named by `CHROMIUM`), Pillow and numpy. The run
  shown is simulated: four barcodes, from the genome of the benchmark
  (`benchmark/README.md`), with different numbers of reads and shares of
  internal priming, so that the samples differ in state:

      benchmark/make-genome.sh genome
      S="python3 -I benchmark/simulate.py genome/genome.fasta"
      $S sim/b1 --seed 1 --shuffle
      $S sim/b2 --seed 2 --reads-per-segment 800 --internal-share 0.2 \
        --joined-share 0.05 --shuffle
      $S sim/b3 --seed 3 --segment-reads 1:60,2:40,3:90 --shuffle
      $S sim/b4 --seed 4 --reads-per-segment 15 --shuffle
      for i in 1 2 3 4; do
        mkdir -p run/fastq_pass/barcode0$i
        cp sim/b$i/reads.fastq.gz run/fastq_pass/barcode0$i/
      done
      flucut -q -r genome/genome.fasta -o out run/fastq_pass
      python3 docs/manual/figures/report_figures.py \
        out/flucut_report.html out/flucut_dashboard.html \
        docs/manual/figures
