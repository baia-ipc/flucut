# (c) 2026 Giorgio Gonnella
"""Common elements of the figures: colours, bars, labels, saving."""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

TEMPLATE = "#3b7fc4"
PRIMER = "#e2a12c"
TAIL = "#8e62c0"
FLANK = "#a3acb6"
BAD = "#d4453c"
OK = "#2c9a5f"
INK = "#1c2530"
MUTED = "#5d6b7a"
LINE = "#d5dbe2"

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "figures")

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "svg.fonttype": "none",
  # (the same file at each drawing)
  "svg.hashsalt": "flucut",
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6,
})


def canvas(width, height, xmax, ymax):
  """A figure without axes, in the units 0..xmax, 0..ymax."""
  fig, ax = plt.subplots(figsize=(width, height))
  ax.set_xlim(0, xmax)
  ax.set_ylim(0, ymax)
  ax.axis("off")
  fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
  return fig, ax


def bar(ax, x0, x1, y, color, h=0.5, text=None, ghost=False, **kw):
  """A stretch of sequence from x0 to x1, centred on y."""
  style = dict(facecolor=color, edgecolor="none", alpha=0.25 if ghost
               else 1.0)
  style.update(kw)
  ax.add_patch(Rectangle((x0, y - h / 2), x1 - x0, h, **style))
  if text:
    ax.text((x0 + x1) / 2, y, text, ha="center", va="center",
            color="white", fontsize=7)


def outline(ax, x0, x1, y, color, h=0.5, **kw):
  """An outlined stretch (a probe, a window)."""
  style = dict(facecolor="none", edgecolor=color, linewidth=1.2)
  style.update(kw)
  ax.add_patch(Rectangle((x0, y - h / 2), x1 - x0, h, **style))


def label(ax, x, y, text, **kw):
  style = dict(ha="left", va="center", fontsize=8, color=INK)
  style.update(kw)
  ax.text(x, y, text, **style)


def bracket(ax, x0, x1, y, text, color=OK, drop=0.18, **kw):
  """A bracket under a stretch, with its meaning."""
  ax.plot([x0, x0, x1, x1], [y + drop, y, y, y + drop], color=color,
          linewidth=1.2, solid_capstyle="butt")
  style = dict(ha="center", va="top", fontsize=7.5, color=color)
  style.update(kw)
  ax.text((x0 + x1) / 2, y - 0.08, text, **style)


def legend(ax, x, y, items, step):
  """Small squares with their meaning, in a row."""
  for i, (color, text) in enumerate(items):
    ax.add_patch(FancyBboxPatch((x + i * step, y - 0.11), 1.6, 0.22,
                                boxstyle="round,pad=0.02",
                                facecolor=color, edgecolor="none"))
    ax.text(x + i * step + 2.4, y, text, va="center", fontsize=7.5,
            color=MUTED)


def save(fig, name):
  """Write the figure as SVG (for the Markdown file) and PDF."""
  os.makedirs(OUT_DIR, exist_ok=True)
  for ext, undated in (("svg", {"Date": None}),
                       ("pdf", {"CreationDate": None})):
    fig.savefig(os.path.join(OUT_DIR, name + "." + ext), metadata=undated)
  plt.close(fig)


def box(ax, x0, x1, y0, y1, text="", color=LINE, fill="white", **kw):
  """A rounded box with a text in it (a step, a group of steps)."""
  ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                              boxstyle="round,pad=0.0,rounding_size=0.12",
                              facecolor=fill, edgecolor=color,
                              linewidth=1.0))
  style = dict(ha="center", va="center", fontsize=7.5, color=INK)
  style.update(kw)
  if text:
    ax.text((x0 + x1) / 2, (y0 + y1) / 2, text, **style)


def plain_axes(ax):
  """Axes of a plot: no top and right border."""
  for side in ("top", "right"):
    ax.spines[side].set_visible(False)
