# flucut brochure

A two-page A4 advertisement of flucut (`brochure.html`), rendered to
`flucut-brochure.pdf`.

- Page 1: the banner, what flucut does, the speed/memory/tolerance
  figures, how it works, and the terminal (banner and colored messages).
- Page 2: the coverage-by-read-origin plot on a real H5N1 Nanopore run,
  a feature comparison with the other trimmers, the live mode, and how
  to get flucut.

## Rebuild

Requires Playwright with Chromium (as the manual figures).

    node pdf.js brochure.html flucut-brochure.pdf

## Images (`img/`)

- `banner.png`: the wide logo (`docs/logo/flucut-banner.png`).
- `terminal_flucut.png`, `read_layout.png`: from the manual figures.
- `coverage.png`: the coverage section of the report of a real Nanopore
  run (SRA `SRR36760543`, segments PB1 and HA, complete H5N1 references
  with the UTRs), rendered as the manual figure of the same name.
- `live_report.png`: the head of the report during a live run.

The numbers come from `benchmark/RESULTS.md` and
`benchmark/COMPARISON.md`.
