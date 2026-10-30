#!/usr/bin/env python3
"""
Diagrams of the libcbaia manual: the components and modules of the
library, and the scheme of multithreaded processing.

Usage: make_diagrams.py [OUTDIR]  (default: the directory of the script)

Writes <name>.pdf (used in the PDF manual) and <name>.png for each
diagram. Requires matplotlib.
"""
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# colors of the components (Okabe-Ito, colorblind-safe), as fill colors
COLORS = {"core": "#E69F00", "ds": "#56B4E9", "io": "#009E73",
          "seq": "#CC79A7", "ext": "#BBBBBB"}
INK = "#222222"
MUTED = "#666666"

def setup(width, height):
  fig, ax = plt.subplots(figsize=(width, height))
  ax.set_axis_off()
  plt.rcParams["font.family"] = "DejaVu Sans"
  return fig, ax

def box(ax, x, y, w, h, text, color, alpha=0.25, bold=False, size=8.5):
  ax.add_patch(FancyBboxPatch((x, y), w, h,
               boxstyle="round,pad=0.02,rounding_size=0.08",
               facecolor=color, alpha=alpha, edgecolor=color, lw=1.2))
  ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
          fontsize=size, color=INK, fontweight="bold" if bold else "normal")

def arrow(ax, x0, y0, x1, y1, color=MUTED):
  ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
              arrowprops=dict(arrowstyle="-|>", color=color, lw=1,
                              shrinkA=2, shrinkB=2))

def save(fig, outdir, name):
  for ext in ("pdf", "png"):
    fig.savefig(os.path.join(outdir, name + "." + ext), bbox_inches="tight",
                dpi=200)
  plt.close(fig)

def components(outdir):
  """the components, their modules, and the main dependencies"""
  fig, ax = setup(9, 6.3)
  W, H, DY = 1.9, 0.42, 0.6
  # component: (x of the frame, columns, modules)
  comps = [
    ("core", 0.0, 2, ["alloc", "logger", "argparse-def", "file",
                       "strutils", "strbuf", "progress", "colors", "test",
                       "numutils"]),
    ("ds", 4.4, 1, ["bitwise", "bloom-filter", "set", "heap",
                     "hyperloglog"]),
    ("io", 6.65, 1, ["file-sorter", "file-merger", "gzip-compressor",
                      "bgzf-compressor", "md5", "html-writer",
                      "html-template", "html-page"]),
  ]
  top = 4.6
  pos = {}
  for name, x, ncols, mods in comps:
    nrows = (len(mods) + ncols - 1) // ncols
    w = ncols * (W + 0.15) + 0.15
    h = nrows * DY + 0.45
    ax.add_patch(FancyBboxPatch((x, top - h), w, h,
                 boxstyle="round,pad=0.02,rounding_size=0.12",
                 facecolor="white", edgecolor=COLORS[name], lw=1.5))
    ax.text(x + 0.12, top - 0.1, name, ha="left", va="top", fontsize=10,
            fontweight="bold", color=INK)
    for k, m in enumerate(mods):
      mx = x + 0.15 + (k % ncols) * (W + 0.15)
      my = top - 0.45 - (k // ncols + 1) * DY + (DY - H)
      box(ax, mx, my, W, H, m, COLORS[name])
      pos[m] = (mx, my)
  # seq: two columns, below
  seq = ["sequtil", "sam-header", "fastq-reader", "bam-reader",
         "fastq-writer", "reads-reader", "ubam-writer", "outfiles",
         "sam-tags", "edit-distance", "minknow-metadata", "whitelist",
         "complexity", "cell-caller"]
  x, ncols = 2.2, 2
  stop = top - 3.55
  nrows = (len(seq) + ncols - 1) // ncols
  w = ncols * (W + 0.15) + 0.15
  h = nrows * DY + 0.45
  ax.add_patch(FancyBboxPatch((x, stop - h), w, h,
               boxstyle="round,pad=0.02,rounding_size=0.12",
               facecolor="white", edgecolor=COLORS["seq"], lw=1.5))
  ax.text(x + 0.12, stop - 0.1, "seq", ha="left", va="top", fontsize=10,
          fontweight="bold", color=INK)
  for k, m in enumerate(seq):
    mx = x + 0.15 + (k % ncols) * (W + 0.15)
    my = stop - 0.45 - (k // ncols + 1) * DY + (DY - H)
    box(ax, mx, my, W, H, m, COLORS["seq"])
    pos[m] = (mx, my)
  def edge(a, b):
    """arrow from module a to module b (a uses b), between box borders"""
    (xa, ya), (xb, yb) = pos[a], pos[b]
    ca = (xa + W / 2, ya + H / 2)
    cb = (xb + W / 2, yb + H / 2)
    if abs(ca[1] - cb[1]) < 1e-6:          # same row
      sx = 1 if cb[0] > ca[0] else -1
      arrow(ax, ca[0] + sx * W / 2, ca[1], cb[0] - sx * W / 2, cb[1])
    else:
      sy = 1 if cb[1] > ca[1] else -1
      arrow(ax, ca[0], ca[1] + sy * H / 2, cb[0], cb[1] - sy * H / 2)
  edge("reads-reader", "bam-reader")
  edge("bam-reader", "fastq-reader")
  edge("cell-caller", "whitelist")
  edge("sam-header", "bam-reader")
  edge("ubam-writer", "sam-tags")
  edge("minknow-metadata", "sam-tags")
  edge("file-merger", "heap")
  edge("html-page", "html-template")
  edge("html-template", "html-writer")
  # external code
  ex, ey = 7.3, stop - h
  box(ax, ex, ey + 1.0, 1.9, 0.42, "ext/zlib", COLORS["ext"], alpha=0.5)
  box(ax, ex, ey + 0.3, 1.9, 0.42, "ext/khash", COLORS["ext"], alpha=0.5)
  ax.text(ex + 0.95, ey + 1.6, "external code", ha="center", fontsize=8.5,
          color=MUTED)
  ax.text(0.0, stop - 0.6, "all modules use\nalloc and logger;\n"
          "zlib: compressors,\nFASTQ and BAM readers;\nkhash: whitelist,\n"
          "cell caller;\noutfiles: BGZF\nend-of-file marker", fontsize=8,
          color=MUTED, va="top")
  ax.set_xlim(-0.1, 9.3)
  ax.set_ylim(stop - h - 0.1, top + 0.1)
  save(fig, outdir, "components")

def threads(outdir):
  """batches read in turn, processed in parallel, written in order"""
  fig, ax = setup(9, 3.0)
  colors = ["#E69F00", "#56B4E9", "#009E73"]
  READ, WRITE = 0.3, 0.35
  process = [1.9, 1.5, 2.2, 1.7, 2.0, 1.6, 1.8, 2.1, 1.5]   # by batch
  n_threads, n_batches = 3, 9
  free = [0.0] * n_threads       # time at which each thread is free
  input_free = 0.0               # time at which the input lock is free
  output_free = 0.0              # time at which the batch can be written
  for b in range(n_batches):
    t = min(range(n_threads), key=lambda k: free[k])
    y = 2.0 - t * 0.75
    r0 = max(free[t], input_free)
    input_free = r0 + READ
    p0 = r0 + READ
    p1 = p0 + process[b]
    w0 = max(p1, output_free)
    output_free = w0 + WRITE
    free[t] = output_free
    box(ax, r0, y, READ, 0.4, "", "#888888", alpha=0.5)
    box(ax, p0, y, process[b], 0.4, f"batch {b}", colors[t], alpha=0.3,
        size=8)
    if w0 > p1 + 0.01:     # waiting for the previous batches
      ax.plot([p1 + 0.03, w0 - 0.03], [y + 0.2, y + 0.2], color=MUTED,
              lw=1, ls=(0, (2, 2)))
    box(ax, w0, y, WRITE, 0.4, "w", colors[t], alpha=0.9, size=8)
  for t in range(n_threads):
    ax.text(-0.15, 2.2 - t * 0.75, f"thread {t + 1}", ha="right",
            va="center", fontsize=9, color=INK)
  # legend
  ly = -0.35
  box(ax, 0.0, ly, 0.3, 0.3, "", "#888888", alpha=0.5)
  ax.text(0.4, ly + 0.15, "read a batch (input lock)", fontsize=8,
          color=MUTED, va="center")
  box(ax, 3.2, ly, 0.6, 0.3, "", colors[0], alpha=0.3)
  ax.text(3.9, ly + 0.15, "process and compress (parallel)", fontsize=8,
          color=MUTED, va="center")
  ax.plot([6.95, 7.35], [ly + 0.15, ly + 0.15], color=MUTED, lw=1,
          ls=(0, (2, 2)))
  ax.text(7.45, ly + 0.15, "wait", fontsize=8, color=MUTED, va="center")
  box(ax, 7.95, ly, 0.35, 0.3, "w", colors[0], alpha=0.9, size=7)
  ax.text(8.4, ly + 0.15, "write in order", fontsize=8, color=MUTED,
          va="center")
  end = max(free)
  ax.annotate("", xy=(end, 2.65), xytext=(0.0, 2.65),
              arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1))
  ax.text(end / 2, 2.72, "time", ha="center", fontsize=8, color=MUTED)
  ax.set_xlim(-1.0, max(end, 9.9) + 0.1)
  ax.set_ylim(-0.5, 2.95)
  save(fig, outdir, "threads")

def main():
  outdir = sys.argv[1] if len(sys.argv) > 1 else \
           os.path.dirname(os.path.abspath(__file__))
  components(outdir)
  threads(outdir)

if __name__ == "__main__":
  main()
