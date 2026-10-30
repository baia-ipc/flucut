# (c) 2026 Giorgio Gonnella
"""Draw the diagrams of the flucut manual.

Writes, in the directory of the script, <name>.pdf (for the PDF manual)
and <name>.png (for the HTML version) of:

  read_layout       the parts of a read of an amplicon and what each
                    trimming policy keeps
  internal_priming  the internal priming site of the HA of KF918483.1
                    (the I348K artefact) and the amplicon it starts

The sequences are those of the primers (src/primers.h) and of the HA
record KF918483.1 (NCBI) at the site (see benchmark/README.md). The
colors are the Okabe-Ito palette (colorblind-safe).

Usage: python3 docs/manual/figures/make_diagrams.py (requires matplotlib)
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
# Okabe-Ito
ORANGE = "#E69F00"
SKY = "#56B4E9"
GREEN = "#009E73"
BLUE = "#0072B2"
VERMILLION = "#D55E00"
PURPLE = "#CC79A7"
GREY = "#999999"


def save(fig, name):
  for ext in ("pdf", "png"):
    fig.savefig(os.path.join(HERE, "%s.%s" % (name, ext)), dpi=200,
                bbox_inches="tight")
  plt.close(fig)


def box(ax, x, w, y, color, label, h=0.6, text_color="white"):
  ax.add_patch(Rectangle((x, y), w, h, facecolor=color, edgecolor="white"))
  if label:
    ax.text(x + w / 2, y + h / 2, label, ha="center", va="center",
            fontsize=8, color=text_color)


def bracket(ax, x0, x1, y, label, color):
  ax.plot([x0, x0, x1, x1], [y - 0.12, y, y, y - 0.12], color=color, lw=1.2)
  ax.text((x0 + x1) / 2, y + 0.08, label, ha="center", va="bottom",
          fontsize=8, color=color)


def read_layout():
  fig, ax = plt.subplots(figsize=(8, 2.6))
  parts = [(0, 1.0, GREY, "flank"),
           (1.0, 1.0, ORANGE, "Uni12\n(terminus)"),
           (2.0, 5.0, BLUE, "template (segment)"),
           (7.0, 1.0, VERMILLION, "Uni13\nterminus"),
           (8.0, 0.8, PURPLE, "tail"),
           (8.8, 1.0, GREY, "flank")]
  for x, w, c, label in parts:
    box(ax, x, w, 0, c, label)
  bracket(ax, 2.0, 7.0, 0.85, "no references (default): kept", BLUE)
  bracket(ax, 1.0, 8.0, 1.45, "with references (terminal): kept",
          VERMILLION)
  ax.text(0, -0.3, "5'", fontsize=8, ha="left", va="top")
  ax.text(9.8, -0.3, "3'", fontsize=8, ha="right", va="top")
  ax.text(4.9, -0.45, "a read of the plus strand of an amplicon (scheme "
          "uni12-mbt13); a read of the minus strand is the reverse "
          "complement", fontsize=7, ha="center", va="top", color="#444444")
  ax.set_xlim(-0.1, 9.9)
  ax.set_ylim(-0.8, 2.0)
  ax.axis("off")
  save(fig, "read_layout")


def internal_priming():
  fig, ax = plt.subplots(figsize=(8, 3.6))
  # the segment and the two amplicons, to scale (1779 bases)
  scale = 9.0 / 1779
  box(ax, 0, 1779 * scale, 2.6, BLUE, "HA segment (KF918483.1, 1779 bases)",
      h=0.45)
  site = 1077 * scale
  ax.plot([site, site], [2.45, 3.25], color=VERMILLION, lw=1.5)
  ax.text(site, 3.3, "Uni12-like site, 1078-1089", fontsize=8,
          ha="center", va="bottom", color=VERMILLION)
  box(ax, 0, 1779 * scale, 1.9, SKY, "amplicon from the terminus",
      h=0.35, text_color="black")
  box(ax, site, (1779 - 1077) * scale, 1.4, ORANGE,
      "amplicon started at the site", h=0.35, text_color="black")
  # the bases at the site
  template = "AGCTATAGCAGG"
  primer = "AGCRAAAGCAGG"
  x0, w = 1.2, 0.45
  ax.text(x0 - 0.2, 0.75, "template", ha="right", va="center", fontsize=8)
  ax.text(x0 - 0.2, 0.25, "primer", ha="right", va="center", fontsize=8)
  for i, (t, p) in enumerate(zip(template, primer)):
    diff = t != p
    box(ax, x0 + i * w, w, 0.55, BLUE, t, h=0.4)
    box(ax, x0 + i * w, w, 0.05, VERMILLION if diff else ORANGE, p, h=0.4)
  # the codons
  for start, label in [(1, "Ala"), (4, "Ile \u2192 Lys"), (7, "Ala")]:
    ax.plot([x0 + start * w + 0.03, x0 + (start + 3) * w - 0.03],
            [-0.1, -0.1], color="#444444", lw=1)
    ax.text(x0 + (start + 1.5) * w, -0.18, label, ha="center", va="top",
            fontsize=7.5, color="#444444")
  ax.text(x0 + 12 * w + 0.3, 0.5,
          "the primer writes its own bases into the reads:\n"
          "a single base changes a codon (I348K)",
          fontsize=8, ha="left", va="center")
  ax.set_xlim(-0.6, 9.4)
  ax.set_ylim(-0.6, 3.7)
  ax.axis("off")
  save(fig, "internal_priming")


if __name__ == "__main__":
  read_layout()
  internal_priming()
