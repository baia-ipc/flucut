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
