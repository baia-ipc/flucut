# (c) 2026 Giorgio Gonnella
"""Figures about the summaries of a sample: internal priming sites,
amplicon lengths. The counts are examples, not data."""

from math import exp
import matplotlib.pyplot as plt
from style import BAD, INK, MUTED, OK, PRIMER, TEMPLATE, plain_axes, save


def sites():
  """Starts of internally primed reads, clustered into sites."""
  counts = {1049: 3, 1051: 14, 1052: 88, 1053: 41, 1055: 6,
            1071: 2, 1073: 9, 1074: 4}
  fig, ax = plt.subplots(figsize=(5.2, 2.3))
  xs = list(range(1044, 1080))
  ax.bar(xs, [counts.get(x, 0) for x in xs], width=0.8, color=PRIMER)
  for first, last, best, total in ((1049, 1055, 1052, 152),
                                   (1071, 1074, 1073, 15)):
    ax.plot([first - 0.5, last + 0.5], [-7, -7], color=INK, linewidth=2,
            clip_on=False)
    ax.text(best, counts[best] + 3, "site %d\n%d reads" % (best, total),
            ha="center", fontsize=7.5, color=INK)
  ax.annotate("more than 5 bases apart:\nanother site", xy=(1063, 2),
              xytext=(1060.5, 45), fontsize=7, color=MUTED, ha="center",
              arrowprops=dict(arrowstyle="-", color=MUTED, linewidth=0.6))
  ax.set_xlabel("position of the first base after the primer")
  ax.set_ylabel("reads that start there")
  ax.set_ylim(0, 115)
  plain_axes(ax)
  fig.tight_layout()
  save(fig, "sites")


def lengths():
  """Histogram of the lengths of complete amplicons of a segment."""
  reference = 1750
  bins = list(range(0, 2000, 10))

  def peak(x, centre, width, height):
    return height * exp(-((x - centre) / width) ** 2)

  values = [peak(b, 1725, 14, 300) + peak(b, 675, 9, 190) +
            peak(b, 980, 60, 12) for b in bins]
  fig, ax = plt.subplots(figsize=(5.2, 2.3))
  limit = 0.9 * reference
  ax.bar([b for b in bins if b < limit],
         [v for b, v in zip(bins, values) if b < limit], width=10,
         align="edge", color=BAD)
  ax.bar([b for b in bins if b >= limit],
         [v for b, v in zip(bins, values) if b >= limit], width=10,
         align="edge", color=TEMPLATE)
  ax.axvline(limit, color=INK, linewidth=0.8, linestyle="--")
  ax.text(limit - 15, 290, "90% of the\nreference", ha="right",
          fontsize=7, va="top")
  ax.text(1725, 310, "length of\nthe reference", ha="center", fontsize=7,
          color=TEMPLATE)
  ax.text(675, 205, "short complete amplicons:\ninternal priming, "
          "deletions", ha="center", fontsize=7, color=BAD)
  ax.set_xlabel("bases between the two primers (bins of 10)")
  ax.set_ylabel("complete amplicons")
  ax.set_xlim(300, 2000)
  ax.set_ylim(0, 360)
  plain_axes(ax)
  fig.tight_layout()
  save(fig, "lengths")
