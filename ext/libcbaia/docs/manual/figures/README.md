# Figures of the libcbaia manual

- `make_diagrams.py`: the diagrams (`components`, `threads`; PDF for the
  manual and PNG); requires matplotlib:

      python3 make_diagrams.py

- `terminal.py` and `screenshot.js`: the screenshots of the terminal
  (`term_*.png`): the commands are run in a pseudo-terminal, so that the
  logger uses colors, and their output is rendered as HTML, then as a
  PNG (requires Playwright with Chromium; `NODE_PATH=$(npm root -g)` for
  a global installation). From `docs/manual/examples`,
  after `make`:

      F=../figures
      python3 $F/terminal.py /tmp/t.html -w 90 -- ./logdemo -p -- ./logdemo -v -n 1000
      node $F/screenshot.js /tmp/t.html $F/term_logger.png
      python3 $F/terminal.py /tmp/t.html -w 90 -- ./logdemo -D -n 1000 -e
      node $F/screenshot.js /tmp/t.html $F/term_debug.png
      python3 $F/terminal.py /tmp/t.html -w 90 -- ./readcount -v data \
        -- ./readcount missing.fastq -- ./readcount -m abc data @1
      node $F/screenshot.js /tmp/t.html $F/term_readcount.png
      python3 $F/terminal.py /tmp/t.html -w 90 -- ./bcdemux -v \
        -w data/barcodes.txt -d 1 -o out data/reads.fastq -- ls out
      node $F/screenshot.js /tmp/t.html $F/term_bcdemux.png
      rm -rf out

- `html_report.png`: the report written by the example program
  `htmlreport`, screenshot of its `<main>` element in a window 900 pixels
  wide (from `docs/manual/examples`, after `make`):

      ./htmlreport /tmp/report.html
      node ../figures/screenshot.js /tmp/report.html ../figures/html_report.png main 900
