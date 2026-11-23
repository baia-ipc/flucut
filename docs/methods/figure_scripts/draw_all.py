#!/usr/bin/env python3
# (c) 2026 Giorgio Gonnella
"""Draw the figures of the methods document into ../figures: SVG, shown
by methods.md, and PDF, for the PDF of the document.

Usage: docs/methods/figure_scripts/draw_all.py   (requires matplotlib)
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import checks  # noqa: E402
import ends  # noqa: E402
import primer_errors  # noqa: E402
import reads  # noqa: E402
import references  # noqa: E402
import run  # noqa: E402
import summaries  # noqa: E402

if __name__ == "__main__":
  reads.amplicons()
  reads.primer_forms()
  reads.split()
  references.probe_location()
  references.junction()
  references.coverage()
  primer_errors.primer_errors()
  run.run_flow()
  run.reference_votes()
  ends.orientations()
  ends.origin()
  ends.adapter_removal()
  ends.from_inside()
  summaries.sites()
  summaries.lengths()
  checks.units()
  checks.trim_verdict()
  checks.simulated_labels()
