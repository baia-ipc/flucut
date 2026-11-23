# (c) 2026 Giorgio Gonnella
"""Figures about the ends of a read: orientation, origin, adapters,
location from inside."""

from style import (BAD, FLANK, INK, LINE, MUTED, OK, PRIMER, TAIL, TEMPLATE,
                   bar, bracket, canvas, label, outline, save)


def orientations():
  """The primers of both ends, searched in both orientations."""
  fig, ax = canvas(6.6, 2.5, 100, 5.0)
  rows = [(3.6, "Read as given", "left primer: not found;\nright primer: "
           "found with 3 errors", False),
          (1.3, "Its reverse complement", "both primers found,\nwithout "
           "errors: this orientation is kept", True)]
  for y, title, verdict, good in rows:
    label(ax, 1, y + 0.75, title, fontweight="bold")
    bar(ax, 3, 67, y, LINE, h=0.5)
    for x0 in (3, 51):
      bar(ax, x0, x0 + 16, y, TEMPLATE, h=0.9, alpha=0.13)
    bar(ax, 3, 9, y, FLANK)
    bar(ax, 61, 67, y, FLANK)
    bar(ax, 13, 57, y, TEMPLATE, text="template")
    label(ax, 70, y, verdict, fontsize=7,
          color="#1d6b41" if good else MUTED)
  # as given: the primers lie as their reverse complements
  bar(ax, 9, 13, 3.6, PRIMER, hatch="////", edgecolor="white")
  bar(ax, 57, 61, 3.6, PRIMER, hatch="////", edgecolor="white")
  bar(ax, 9, 13, 1.3, PRIMER)
  bar(ax, 57, 61, 1.3, PRIMER)
  label(ax, 11, 0.72, "left\nprimer", fontsize=6.5, ha="center", va="top")
  label(ax, 59, 0.72, "right\nprimer", fontsize=6.5, ha="center",
        va="top")
  bracket(ax, 3, 19, 2.85, "window W: the left\nprimer begins here",
          color=TEMPLATE, fontsize=6.5)
  bracket(ax, 51, 67, 2.85, "window W: the right\nprimer ends here",
          color=TEMPLATE, fontsize=6.5)
  save(fig, "orientations")


def origin():
  """The origin of an end, from where its probe lies."""
  fig, ax = canvas(6.6, 2.3, 100, 4.6)
  label(ax, 1, 4.3, "Reference segment", fontweight="bold")
  bar(ax, 5, 95, 3.6, TEMPLATE, h=0.3, alpha=0.45)
  bar(ax, 5, 13, 3.6, OK, h=0.6, alpha=0.3)
  bar(ax, 87, 95, 3.6, OK, h=0.6, alpha=0.3)
  label(ax, 86, 4.25, "terminal windows (60 bases), at both ends",
        fontsize=6.5, color="#1d6b41", ha="right")
  cases = [(6, "TERMINAL", "#1d6b41",
            "the probe lies within 60 bases\nof the terminus"),
           (44, "INTERNAL", BAD,
            "the probe lies further inside:\ninternal priming"),
           (None, "UNKNOWN", MUTED, "the probe is not located")]
  x_text = [5, 40, 76]
  for (x, name, color, text), xt in zip(cases, x_text):
    bar(ax, xt, xt + 3, 1.6, PRIMER)
    outline(ax, xt + 3, xt + 11, 1.6, INK, h=0.6)
    bar(ax, xt + 3, xt + 22, 1.6, TEMPLATE, alpha=0.8)
    outline(ax, xt + 3, xt + 11, 1.6, INK, h=0.6)
    label(ax, xt, 0.85, name, fontsize=8, color=color,
          fontweight="bold")
    label(ax, xt, 0.35, text, fontsize=6.5, color=MUTED)
    if x is not None:
      outline(ax, x, x + 8, 3.6, INK, h=0.5)
      ax.annotate("", xy=(x + 4, 3.3), xytext=(xt + 7, 1.95),
                  arrowprops=dict(arrowstyle="->", color=INK,
                                  linewidth=0.8))
  label(ax, 16.5, 2.25, "probe: the 32 bases after the primer",
        fontsize=6.5, color=INK)
  save(fig, "origin")


def adapter_removal():
  """An end without a primer: adapter removed, second analysis."""
  fig, ax = canvas(6.6, 2.4, 100, 4.8)
  rows = [(4.0, "First analysis", "no primer is found at this end:\ntoo "
           "little of it is left"),
          (2.4, "Adapter", "an adapter is found within the\nwindow and "
           "removed"),
          (0.8, "Second analysis", "a primer is accepted within 5 "
           "bases\nof the adapter, also overhanging")]
  for y, title, text in rows:
    label(ax, 1, y + 0.45, title, fontweight="bold")
    label(ax, 62, y, text, fontsize=7, color=MUTED)
  bar(ax, 3, 21, 4.0, FLANK, text="adapter + barcode")
  bar(ax, 21, 23.5, 4.0, PRIMER)
  bar(ax, 23.5, 58, 4.0, TEMPLATE, text="template")
  bar(ax, 3, 21, 2.4, FLANK, ghost=True)
  outline(ax, 3, 21, 2.4, BAD, h=0.5)
  bar(ax, 21, 23.5, 2.4, PRIMER)
  bar(ax, 23.5, 58, 2.4, TEMPLATE)
  bar(ax, 15, 21, 0.8, PRIMER, ghost=True)
  bar(ax, 21, 23.5, 0.8, PRIMER)
  outline(ax, 15, 23.5, 0.8, BAD, h=0.5)
  bar(ax, 23.5, 58, 0.8, TEMPLATE)
  bracket(ax, 23.5, 58, 0.38, "kept", color=OK)
  save(fig, "adapter_removal")


def from_inside():
  """An end located by probes from further inside the read."""
  fig, ax = canvas(6.6, 3.0, 100, 6.0)
  label(ax, 1, 5.7, "The probe next to the primer is not located: probes "
        "further inside are tried", fontweight="bold")
  bar(ax, 5, 9, 4.9, PRIMER)
  bar(ax, 9, 70, 4.9, TEMPLATE)
  for i in range(3):
    x = 9 + 9 * i
    outline(ax, x, x + 9, 4.9, INK if i == 2 else MUTED, h=0.62,
            linestyle="-" if i == 2 else "--")
    label(ax, x + 4.5, 4.9, str(i + 1), fontsize=7, ha="center",
          color="white")
  label(ax, 37, 4.9, "probe 3 is located", fontsize=7, color="white")
  bracket(ax, 9, 27, 4.4, "bases skipped", color=MUTED)

  label(ax, 1, 3.3, "The bases skipped continue the reference",
        fontweight="bold")
  bar(ax, 5, 95, 2.6, TEMPLATE, h=0.25, alpha=0.45)
  outline(ax, 48, 57, 2.6, INK, h=0.5)
  bar(ax, 30, 48, 2.6, TEMPLATE, h=0.25)
  ax.plot([30], [2.6], marker="v", color=BAD, markersize=6)
  label(ax, 30, 3.0, "position of the end: that of the probe, moved out "
        "by the bases skipped", fontsize=6.5, color=BAD, ha="left")

  label(ax, 1, 1.6, "They do not: a short flank at the terminus, then a "
        "deletion", fontweight="bold")
  bar(ax, 5, 95, 0.9, TEMPLATE, h=0.25, alpha=0.45)
  bar(ax, 5, 9, 0.9, TEMPLATE, h=0.25)
  outline(ax, 48, 57, 0.9, INK, h=0.5)
  ax.plot([5], [0.9], marker="v", color=BAD, markersize=6)
  bar(ax, 9, 48, 0.9, LINE, h=0.5, alpha=0.7)
  label(ax, 28.5, 0.9, "deleted in the read", fontsize=6.5, ha="center",
        color=MUTED)
  label(ax, 5, 0.3, "the first bases skipped fit the terminus: the end "
        "is there", fontsize=6.5, color=BAD)
  save(fig, "from_inside")
