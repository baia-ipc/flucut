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
