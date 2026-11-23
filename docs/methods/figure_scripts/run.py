# (c) 2026 Giorgio Gonnella
"""Figures about a whole run: its flow, the choice of the references."""

import matplotlib.pyplot as plt
from style import (INK, MUTED, OK, PRIMER, TAIL, TEMPLATE, box, canvas,
                   label, plain_axes, save)

SOFT = "#F3F6FB"


def run_flow():
  """The steps of a run, nested as its loops are."""
  fig, ax = canvas(6.6, 3.5, 100, 9.2)
  box(ax, 1, 99, 0.2, 9.0, color=MUTED)
  label(ax, 2.5, 8.6, "Run", fontweight="bold")
  box(ax, 14, 55, 8.25, 8.9, "k-mer index of the references (D1)",
      fill=SOFT)
  box(ax, 57, 97, 8.25, 8.9, "error percentage allowed in a primer (A2)",
      fill=SOFT)
  box(ax, 3, 97, 0.5, 7.9, color=MUTED)
  label(ax, 4.5, 7.5, "Each readset (a sample)", fontweight="bold")
  box(ax, 40, 95, 7.15, 7.8, "closest reference of each segment (A3)",
      fill=SOFT)
  box(ax, 5, 95, 1.7, 6.8, color=MUTED)
  label(ax, 6.5, 6.4, "Each read", fontweight="bold")
  box(ax, 30, 93, 6.05, 6.7, "split, if it holds several amplicons (A4)",
      fill=SOFT)
  box(ax, 7, 93, 1.95, 5.7, color=MUTED)
  label(ax, 8.5, 5.3, "Each part", fontweight="bold")
  steps = [("primers of\nboth ends\n(A5, A6, A7)", TEMPLATE),
           ("origin of\neach end\n(A9, A10)", TEMPLATE),
           ("part to\nkeep\n(A11)", TEMPLATE),
           ("adapters,\nsecond analysis\n(A12)", PRIMER),
           ("positions from\ninside, junction\n(A13, A14)", TAIL),
           ("record\nwritten,\ncounts", OK)]
  x = 9.0
  for i, (text, color) in enumerate(steps):
    box(ax, x, x + 12.6, 2.35, 4.75, text, color=color, fontsize=7)
    if i < len(steps) - 1:
      ax.annotate("", xy=(x + 14.0, 3.55), xytext=(x + 12.7, 3.55),
                  arrowprops=dict(arrowstyle="->", color=INK,
                                  linewidth=0.8))
    x += 14.0
  box(ax, 30, 95, 0.7, 1.5, "internal priming sites (A15) · amplicon "
      "lengths (A16)\ncoverage profiles (A17) · reports (D3)",
      fill=SOFT, fontsize=7)
  save(fig, "run_flow")


def reference_votes():
  """Votes of the probes for the references of five segments."""
  fig, ax = plt.subplots(figsize=(6.2, 2.4))
  # shares of the votes of a segment: subtype A, B, C (an example)
  rows = [("HA", (0.97, 0.02, 0.01), "clear: reference A", 900),
          ("NA", (0.88, 0.09, 0.03), "clear: reference A", 700),
          ("NS", (0.38, 0.52, 0.10),
           "no clear lead, A is a close contender: follows A", 400),
          ("NP", (0.15, 0.50, 0.35),
           "no clear lead, A is not close: most votes, B", 500),
          ("PA", (0.5, 0.3, 0.2), "fewer than 20 votes: follows A", 12)]
  colors = (TEMPLATE, PRIMER, TAIL)
  for i, (segment, shares, verdict, votes) in enumerate(rows):
    y = len(rows) - 1 - i
    left = 0
    for share, color in zip(shares, colors):
      ax.barh(y, 100 * share, left=left, height=0.55, color=color)
      left += 100 * share
    ax.text(104, y, verdict, va="center", fontsize=7.5, color=INK)
    ax.text(-3, y - 0.36, "%d votes" % votes, ha="right", fontsize=6,
            color=MUTED)
  ax.axvline(60, color=INK, linewidth=0.8, linestyle="--")
  ax.text(60, len(rows) - 0.45, "60%", ha="center", fontsize=7,
          color=INK)
  ax.set_yticks(range(len(rows)))
  ax.set_yticklabels([r[0] for r in reversed(rows)])
  ax.set_xlim(0, 100)
  ax.set_xticks([0, 20, 40, 60, 80, 100])
  ax.set_xlabel("share of the votes of the segment (%)")
  for i, (name, color) in enumerate(zip("ABC", colors)):
    ax.text(104 + 30 * i, -1.45, "■ subtype " + name, fontsize=7,
            color=color, clip_on=False)
  ax.text(104, -0.85, "references of the segment:", fontsize=7,
          color=MUTED, clip_on=False)
  plain_axes(ax)
  fig.subplots_adjust(left=0.1, right=0.52, top=0.93, bottom=0.2)
  save(fig, "reference_votes")
