# (c) 2026 Giorgio Gonnella
"""Figures about the checks: amplicons built from alignments, verdicts
on the trimming, labels of a simulated read."""

from style import (BAD, FLANK, INK, LINE, MUTED, OK, PRIMER, TAIL, TEMPLATE,
                   bar, bracket, canvas, label, legend, save)


def link(ax, r0, r1, g0, g1, y_read, y_ref):
  """An alignment: a stretch of the read on a stretch of the reference."""
  ax.fill([r0, r1, g1, g0], [y_read, y_read, y_ref, y_ref],
          color=TEMPLATE, alpha=0.18, linewidth=0)


def units():
  """Alignments of a read joined into amplicons, or not."""
  fig, ax = canvas(6.6, 2.9, 100, 5.8)
  label(ax, 1, 5.5, "One amplicon with a deletion: two alignments, no "
        "primer between them", fontweight="bold")
  bar(ax, 20, 24, 4.8, PRIMER)
  bar(ax, 24, 56, 4.8, TEMPLATE)
  bar(ax, 56, 60, 4.8, PRIMER)
  ax.plot([40, 40], [4.5, 5.1], color="white", linewidth=1.2)
  label(ax, 62, 4.8, "read", fontsize=7, color=MUTED)
  bar(ax, 5, 95, 3.5, TEMPLATE, h=0.25, alpha=0.45)
  link(ax, 24, 40, 5, 21, 4.55, 3.63)
  link(ax, 40, 56, 79, 95, 4.55, 3.63)
  bar(ax, 21, 79, 3.5, LINE, h=0.5, alpha=0.8)
  label(ax, 50, 3.5, "deletion: at least 100 bases of the reference",
        fontsize=6.5, ha="center", color=MUTED)

  label(ax, 1, 2.6, "Two amplicons joined: a primer lies between the two "
        "alignments", fontweight="bold")
  bar(ax, 10, 14, 1.9, PRIMER)
  bar(ax, 14, 44, 1.9, TEMPLATE)
  bar(ax, 44, 48, 1.9, PRIMER)
  bar(ax, 48, 52, 1.9, PRIMER)
  bar(ax, 52, 82, 1.9, TEMPLATE)
  bar(ax, 82, 86, 1.9, PRIMER)
  label(ax, 88, 1.9, "read", fontsize=7, color=MUTED)
  bar(ax, 5, 45, 0.6, TEMPLATE, h=0.25, alpha=0.45)
  bar(ax, 55, 95, 0.6, TEMPLATE, h=0.25, alpha=0.45)
  link(ax, 14, 44, 10, 40, 1.65, 0.73)
  link(ax, 52, 82, 58, 88, 1.65, 0.73)
  label(ax, 25, 0.2, "a segment", fontsize=6.5, ha="center", color=MUTED)
  label(ax, 75, 0.2, "the same or another segment", fontsize=6.5,
        ha="center", color=MUTED)
  save(fig, "units")


def trim_verdict():
  """The verdict on the bases removed at an end, by their number."""
  fig, ax = canvas(6.6, 2.2, 100, 4.4)
  tau = 12   # length of the terminus, as an example

  def x_of(delta):
    return 58 + 2.4 * delta

  rows = [(3.0, "End with an internal origin",
           [(-22, -2.5, BAD, "primer bases left"),
            (-2.5, 2.5, OK, "whole primer\nremoved"),
            (2.5, 16, PRIMER, "template removed")]),
          (1.1, "End with a terminal origin (terminus of 12 bases)",
           [(-22, -tau - 2.5, BAD, "tail or adapter\nbases left"),
            (-tau - 2.5, -tau + 2.5, OK, "terminus\nkept"),
            (-tau + 2.5, -2.5, TAIL, "part of the terminus removed"),
            (-2.5, 2.5, OK, "whole primer\nremoved"),
            (2.5, 16, PRIMER, "template removed")])]
  for y, title, zones in rows:
    label(ax, 1, y + 0.62, title, fontweight="bold")
    for d0, d1, color, text in zones:
      bar(ax, x_of(d0), x_of(d1), y, color, h=0.55, alpha=0.85)
      label(ax, (x_of(d0) + x_of(d1)) / 2, y, text, fontsize=6,
            ha="center", color="white")
  for d in (-12, -2, 0, 2):
    ax.plot([x_of(d), x_of(d)], [0.62, 0.8], color=INK, linewidth=0.7)
    label(ax, x_of(d), 0.47, ("%+d" % d) if d else "0", fontsize=6.5,
          ha="center")
  label(ax, 58, 0.02, "Δ: bases removed at the read end, less the bases "
        "outside the template there", fontsize=7, ha="center",
        color=MUTED, va="bottom")
  save(fig, "trim_verdict")


def simulated_labels():
  """The label of each base of a simulated read."""
  fig, ax = canvas(6.6, 1.9, 100, 3.8)
  rows = [(3.0, "Amplicon primed at the termini", PRIMER, "p"),
          (1.6, "Amplicon primed at the internal site of HA", BAD, "i")]
  for y, title, color, letter in rows:
    label(ax, 1, y + 0.5, title, fontweight="bold")
    start = 3 if letter == "p" else 33
    bar(ax, start, start + 8, y, FLANK, text="f")
    bar(ax, start + 8, start + 13, y, color, text=letter)
    bar(ax, start + 13, 80, y, TEMPLATE, text="T")
    bar(ax, 80, 85, y, PRIMER, text="p")
    bar(ax, 85, 90, y, TAIL, text="t")
    bar(ax, 90, 98, y, FLANK, text="f")
  legend(ax, 2, 0.75, [(FLANK, "f  flank"), (TAIL, "t  primer tail"),
                       (TEMPLATE, "T  template")], 30)
  legend(ax, 2, 0.25, [(PRIMER, "p  primer at a terminus"),
                       (BAD, "i  primer at the internal site")], 30)
  save(fig, "simulated_labels")
