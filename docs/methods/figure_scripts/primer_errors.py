# (c) 2026 Giorgio Gonnella
"""Plot: the share of primers with more errors than allowed, by the
error rate of the reads (binomial model, as in section A2)."""

from math import comb
import matplotlib.pyplot as plt
from style import BAD, INK, MUTED, TEMPLATE, save


def share_above(length, allowed, rate):
  """Probability of more than <allowed> errors in <length> bases."""
  within = sum(comb(length, k) * rate ** k * (1 - rate) ** (length - k)
               for k in range(allowed + 1))
  return 1 - within


def primer_errors():
  rates = [i / 1000 for i in range(0, 151)]
  fig, ax = plt.subplots(figsize=(4.6, 2.7))
  curves = [(2, "limit 20%: more than 2 errors", BAD),
            (3, "limit 30%: more than 3 errors", TEMPLATE)]
  for allowed, text, color in curves:
    ax.plot([100 * r for r in rates],
            [100 * share_above(12, allowed, r) for r in rates],
            color=color, linewidth=1.8, label=text)
    y = 100 * share_above(12, allowed, 0.08)
    ax.plot([8], [y], marker="o", color=color, markersize=4)
    ax.annotate("%.1f%%" % y, xy=(8, y), xytext=(6.2, y + 1.6),
                fontsize=7.5, color=color)
  ax.axvline(8, color=MUTED, linewidth=0.6, linestyle=":")
  ax.set_xlabel("errors per base in the reads (%)")
  ax.set_ylabel("Uni12 primers missed (%)")
  ax.set_xlim(0, 15)
  ax.set_ylim(0, None)
  ax.legend(frameon=False, fontsize=7.5, loc="upper left",
            labelcolor=INK)
  for side in ("top", "right"):
    ax.spines[side].set_visible(False)
  fig.tight_layout()
  save(fig, "primer_errors")
