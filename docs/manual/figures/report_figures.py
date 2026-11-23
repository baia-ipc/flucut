#!/usr/bin/env python3
"""
Figures of the HTML report and of the dashboard, for the manual.

Usage:
  report_figures.py <report.html> <dashboard.html> <outdir>

The report is rendered by a headless Chromium (the program named by the
environment variable CHROMIUM, default: chromium) as one tall image,
which is cut at the lines above its headings and at the gaps between
its cards: report_overview.png (head and run overview), report_ha.png
(the HA in all samples), report_sample.png, report_segments.png and
report_reads.png (the first sample). The top of the dashboard is saved
as dashboard.png. Requires Pillow and numpy.
"""
import os
import subprocess
import sys
import tempfile
import numpy as np
from PIL import Image

WIDTH = 1200
ACCENT = (248, 74, 42)  # the line above a heading of level 2


def render(html, width, height):
  """A screenshot of the top <height> pixels of the page <html>."""
  with tempfile.TemporaryDirectory() as tmp:
    png = os.path.join(tmp, "page.png")
    subprocess.run(
      [os.environ.get("CHROMIUM", "chromium"), "--headless", "--no-sandbox",
       "--disable-gpu", "--hide-scrollbars",
       f"--window-size={width},{height}", f"--screenshot={png}",
       "file://" + os.path.abspath(html)],
      check=True, stderr=subprocess.DEVNULL)
    return Image.open(png).convert("RGB")


def heading_lines(pixels):
  """The rows of the lines above the headings of level 2."""
  accent = (pixels == ACCENT).all(axis=2).sum(axis=1)
  rows = [y for y in range(len(accent)) if accent[y] > WIDTH // 2]
  return [y for i, y in enumerate(rows) if i == 0 or y > rows[i - 1] + 1]


def gaps(pixels, start, end, min_height=10):
  """The white bands (first row, row after the last) from <start> to <end>."""
  white = (pixels[start:end] == 255).all(axis=(1, 2))
  bands, first = [], None
  for y, is_white in enumerate(white):
    if is_white and first is None:
      first = y
    if not is_white and first is not None:
      if y - first >= min_height:
        bands.append((start + first, start + y))
      first = None
  if first is not None:
    bands.append((start + first, end))
  return bands


def main(report, dashboard, outdir):
  page = render(report, WIDTH, 7000)
  pixels = np.array(page)
  lines = heading_lines(pixels)
  # the head, then the run overview, the HA and the first sample
  overview, ha, sample, following = lines[0], lines[1], lines[2], lines[3]

  def save(name, top, bottom):
    page.crop((0, top, WIDTH, bottom)).save(
      os.path.join(outdir, name + ".png"), optimize=True)

  save("report_overview", 0, ha - 10)
  save("report_ha", ha + 2, ha + 600)
  # the cards of the sample: reads and sites, segments, read lengths
  bands = gaps(pixels, sample + 150, following)
  cards = [(a[1], b[0]) for a, b in zip(bands, bands[1:])]
  save("report_sample", sample + 2, cards[0][1] + 4)
  save("report_segments", cards[1][0] - 4, cards[1][0] + 1020)
  save("report_reads", cards[2][0] - 4, cards[2][1] + 4)
  top = render(dashboard, 1300, 900)
  rows = np.where(~(np.array(top) == 255).all(axis=(1, 2)))[0]
  top.crop((0, 0, 1300, min(int(rows[-1]) + 12, 900))).save(
    os.path.join(outdir, "dashboard.png"), optimize=True)


if __name__ == "__main__":
  if len(sys.argv) != 4:
    sys.exit(__doc__)
  main(*sys.argv[1:])
