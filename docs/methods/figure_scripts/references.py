# (c) 2026 Giorgio Gonnella
"""Figures about the references: location of a probe, junction of an
internal deletion, coverage by read origin."""

import matplotlib.pyplot as plt
from style import (BAD, INK, LINE, MUTED, OK, PRIMER, TEMPLATE, bar,
                   bracket, canvas, label, outline, save)


def probe_location():
  """A probe located through its seeds."""
  fig, ax = canvas(6.6, 2.3, 100, 4.6)
  label(ax, 1, 4.3, "Probe (32 bases) and its 10-mers", fontweight="bold")
  bar(ax, 20, 44, 3.7, TEMPLATE, text="probe")
  ax.plot([33], [3.7], marker="x", color=BAD, markersize=7,
          markeredgewidth=1.8)
  label(ax, 46, 3.7, "a sequencing error: the 10-mers over it\nare not "
        "in the references", fontsize=7, color=MUTED)
  seeds = [20, 36.5]
  for x in seeds:
    bar(ax, x, x + 7.5, 3.2, PRIMER, h=0.12)
    ax.fill([x, x + 7.5, x + 42.5, x + 35], [3.14, 3.14, 1.35, 1.35],
            color=PRIMER, alpha=0.18, linewidth=0)
  label(ax, 1, 1.75, "Reference", fontweight="bold")
  bar(ax, 5, 97, 1.2, TEMPLATE, h=0.3, alpha=0.45)
  for x in seeds:
    bar(ax, x + 35, x + 42.5, 1.2, PRIMER, h=0.3)
  outline(ax, 51, 83, 1.2, INK, h=0.62)
  bracket(ax, 51, 83, 0.72, "window aligned: around the diagonal of the "
          "seeds", color=INK)
  label(ax, 66, 2.45, "10-mers found in the reference (seeds),\non the "
        "same diagonal", fontsize=7, color="#b07a12")
  save(fig, "probe_location")


def junction():
  """The bisection for the junction of an internal deletion."""
  fig, ax = plt.subplots(figsize=(5.2, 3.3))
  ax.axvspan(420, 1780, color=LINE, alpha=0.6, linewidth=0)
  ax.plot([0, 420], [0, 420], color=TEMPLATE, linewidth=2.2)
  ax.plot([1780, 2300], [420, 940], color=TEMPLATE, linewidth=2.2)
  ax.text(1100, 150, "deleted in the read", ha="center", color=MUTED,
          fontsize=8)
  ax.text(250, 180, "left flank", color=TEMPLATE, fontsize=8)
  ax.text(2080, 800, "right flank", color=TEMPLATE, fontsize=8,
          ha="right")
  # the diagonals continued: where a probe would lie on each flank
  ax.plot([420, 940], [420, 940], color=TEMPLATE, linewidth=0.8,
          linestyle=":")
  ax.plot([1360, 1780], [0, 420], color=TEMPLATE, linewidth=0.8,
          linestyle=":")
  probes = [(470, "right", "1"), (235, "left", "2"), (352, "left", "3")]
  for y, flank, n in probes:
    x = y if flank == "left" else y + 1360
    ax.plot([x], [y], marker="o", color=PRIMER, markersize=9,
            markeredgecolor=INK, markeredgewidth=0.6)
    ax.text(x, y, n, ha="center", va="center", fontsize=6.5)
  ax.axhspan(352, 502, xmin=0, xmax=1, color=OK, alpha=0.18, linewidth=0)
  ax.text(1100, 427, "aligned to both flanks: the junction is where\n"
          "the two alignments have the fewest errors", color="#1d6b41",
          fontsize=7, ha="center", va="center")
  ax.text(1100, 690, "probes 1, 2, 3: each lies on one flank and\n"
          "halves the part not yet assigned", color=INK, fontsize=7,
          ha="center")
  ax.set_xlim(0, 2300)
  ax.set_ylim(0, 940)
  ax.set_xlabel("position in the reference")
  ax.set_ylabel("position in the template of the read")
  for side in ("top", "right"):
    ax.spines[side].set_visible(False)
  fig.tight_layout()
  save(fig, "junction")


def coverage():
  """Spans of reads by origin, and the depth they give."""
  fig, (top, low) = plt.subplots(2, 1, figsize=(5.6, 3.1), sharex=True,
                                 gridspec_kw=dict(height_ratios=[1.25, 1]))
  site = 1050
  for i in range(4):
    top.barh(7 - i, 1700, left=0, height=0.6, color=TEMPLATE)
  for i in range(3):
    top.barh(3 - i, 1700 - site, left=site, height=0.6, color=PRIMER)
  top.text(20, 7.75, "reads that start at the terminus (terminal)",
           fontsize=7.5, color=TEMPLATE)
  top.text(site - 25, 2, "reads that start at an\ninternal priming site",
           fontsize=7.5, color="#b07a12", ha="right", va="center")
  top.set_ylim(0.3, 8.4)
  top.set_yticks([])
  top.set_ylabel("spans")
  xs = [0, site, site, 1700]
  low.fill_between(xs, [4, 4, 4, 4], color=TEMPLATE, linewidth=0)
  low.fill_between(xs, [4, 4, 7, 7], [4, 4, 4, 4], color=PRIMER,
                   linewidth=0)
  low.annotate("step at the site", xy=(site, 5.5), xytext=(site - 420, 6),
               fontsize=7.5, color=INK,
               arrowprops=dict(arrowstyle="->", color=INK, linewidth=0.8))
  low.set_ylim(0, 8)
  low.set_xlim(0, 1700)
  low.set_ylabel("depth")
  low.set_xlabel("position in the reference")
  for ax in (top, low):
    for side in ("top", "right"):
      ax.spines[side].set_visible(False)
  fig.tight_layout()
  save(fig, "coverage")
