# (c) 2026 Giorgio Gonnella
"""Figures about a read: its amplicon, its primers, its splitting."""

from style import (BAD, FLANK, INK, MUTED, OK, PRIMER, TAIL, TEMPLATE, bar,
                   bracket, canvas, label, legend, outline, save)


def amplicons():
  """An amplicon from the termini and an internally primed one."""
  fig, ax = canvas(6.6, 2.75, 100, 5.6)
  label(ax, 2, 5.25, "Segment", fontweight="bold")
  bar(ax, 10, 90, 4.7, TEMPLATE, h=0.32, alpha=0.45)
  bar(ax, 10, 14, 4.7, INK, h=0.32)
  bar(ax, 86, 90, 4.7, INK, h=0.32)
  bar(ax, 42, 46, 4.7, BAD, h=0.32, alpha=0.75)
  label(ax, 10, 4.45, "conserved terminus\n(12 bases)", fontsize=7,
        color=MUTED, va="top")
  label(ax, 90, 4.45, "conserved terminus\n(13 bases)", fontsize=7,
        color=MUTED, ha="right", va="top")
  label(ax, 44, 4.45, "site that resembles\na terminus", fontsize=7,
        color=BAD, ha="center", va="top")

  label(ax, 2, 3.65, "Amplicon primed at the termini", fontweight="bold")
  bar(ax, 5, 10, 3.1, TAIL)
  bar(ax, 10, 14, 3.1, PRIMER)
  bar(ax, 14, 86, 3.1, TEMPLATE, text="template")
  bar(ax, 86, 90, 3.1, PRIMER)
  bar(ax, 90, 95, 3.1, TAIL)
  bracket(ax, 10, 90, 2.62,
          "kept: the termini are those of the segment, only the tails "
          "are removed")

  label(ax, 2, 1.75, "Internally primed amplicon", fontweight="bold")
  bar(ax, 37, 42, 1.2, TAIL)
  bar(ax, 42, 46, 1.2, PRIMER, edgecolor=BAD, linewidth=1.6)
  bar(ax, 46, 86, 1.2, TEMPLATE, text="template")
  bar(ax, 86, 90, 1.2, PRIMER)
  bar(ax, 90, 95, 1.2, TAIL)
  label(ax, 36, 1.2, "bases of the primer,\nnot of the template",
        fontsize=7, color=BAD, ha="right")
  bracket(ax, 46, 90, 0.72, "kept: the whole primer is removed at the "
          "internal start")
  legend(ax, 2, 0.12, [(TEMPLATE, "template"), (PRIMER, "primer terminus"),
                       (TAIL, "primer tail")], 22)
  save(fig, "amplicons")


def primer_forms():
  """The three forms of a primer at the start of a read."""
  fig, ax = canvas(6.6, 2.1, 100, 4.2)
  rows = [(3.5, "Whole primer after a flank",
           "edit distance: substitutions,\ninsertions, deletions"),
          (2.1, "Truncated primer after a flank",
           "at least half of it and 8 bases;\none error fewer allowed"),
          (0.7, "Read that begins inside the primer",
           "tried only if no other form is found;\nerrors in proportion "
           "to the part in the read")]
  for y, title, note in rows:
    label(ax, 1, y + 0.42, title, fontweight="bold")
    label(ax, 64, y, note, fontsize=7, color=MUTED)
  bar(ax, 3, 14, 3.5, FLANK)
  bar(ax, 14, 26, 3.5, PRIMER, text="primer")
  bar(ax, 26, 60, 3.5, TEMPLATE, text="template")
  bar(ax, 3, 14, 2.1, FLANK)
  bar(ax, 14, 19, 2.1, PRIMER, ghost=True)
  bar(ax, 19, 26, 2.1, PRIMER)
  bar(ax, 26, 60, 2.1, TEMPLATE)
  label(ax, 16.5, 1.72, "missing", fontsize=6.5, color=MUTED, ha="center")
  bar(ax, 9, 14, 0.7, PRIMER, ghost=True)
  bar(ax, 14, 21, 0.7, PRIMER)
  bar(ax, 21, 60, 0.7, TEMPLATE)
  ax.plot([14, 14], [0.3, 1.1], color=INK, linewidth=1, linestyle=":")
  label(ax, 14, 0.18, "start of the read", fontsize=6.5, color=MUTED,
        ha="center")
  save(fig, "primer_forms")


def split():
  """The two kinds of junction at which a read is split."""
  fig, ax = canvas(6.6, 3.6, 100, 7.0)
  label(ax, 1, 6.7, "A pair: the end of an amplicon, then the start of "
        "another", fontweight="bold")
  bar(ax, 3, 7, 6.0, PRIMER)
  bar(ax, 7, 40, 6.0, TEMPLATE, text="amplicon 1")
  bar(ax, 40, 44, 6.0, PRIMER)
  bar(ax, 44, 52, 6.0, FLANK)
  bar(ax, 52, 56, 6.0, PRIMER)
  bar(ax, 56, 93, 6.0, TEMPLATE, text="amplicon 2")
  bar(ax, 93, 97, 6.0, PRIMER)
  label(ax, 42, 6.42, "END", fontsize=7, ha="center", color=MUTED)
  label(ax, 54, 6.42, "START", fontsize=7, ha="center", color=MUTED)
  bracket(ax, 44, 52, 5.55, "at most 100 bases:\ndropped", color=MUTED)
  for x in (44, 52):
    ax.plot([x, x], [5.6, 6.4], color=BAD, linewidth=1.4)
  label(ax, 60, 5.3, "at most 2 errors in the two matches together;\n"
        "not split if the references show that the read\nruns through "
        "the pair", fontsize=7, color=MUTED)

  label(ax, 1, 4.3, "A lone match, with references: does the read stop "
        "following the reference there?", fontweight="bold")
  bar(ax, 3, 7, 3.4, PRIMER)
  bar(ax, 7, 46, 3.4, TEMPLATE)
  bar(ax, 46, 50, 3.4, TAIL)
  bar(ax, 50, 54, 3.4, PRIMER)
  bar(ax, 54, 93, 3.4, TEMPLATE)
  bar(ax, 93, 97, 3.4, PRIMER)
  label(ax, 14, 3.4, "amplicon 1", fontsize=7, ha="center", color="white")
  label(ax, 80, 3.4, "amplicon 2", fontsize=7, ha="center", color="white")
  label(ax, 52, 3.86, "START", fontsize=7, ha="center", color=MUTED)
  outline(ax, 54, 62, 3.4, INK, h=0.62)
  label(ax, 58, 3.86, "inner probe", fontsize=7, ha="left")
  for x0 in (38, 30, 22):
    outline(ax, x0, x0 + 8, 3.4, INK, h=0.62, linestyle="--")
  label(ax, 34, 3.86, "outer probes, at three distances", fontsize=7,
        ha="center")
  ax.plot([46, 46], [2.95, 3.85], color=BAD, linewidth=1.4)
  label(ax, 46, 2.85, "cut outside the whole primer", fontsize=7,
        color=BAD, ha="center", va="top")
  # where the probes lie in the references
  label(ax, 1, 1.7, "segment of\namplicon 2", fontsize=7, color=MUTED)
  bar(ax, 16, 97, 1.7, TEMPLATE, h=0.2, alpha=0.45)
  outline(ax, 22, 30, 1.7, INK, h=0.4)
  label(ax, 22, 2.12, "the inner probe is located here", fontsize=7)
  label(ax, 1, 0.7, "segment of\namplicon 1", fontsize=7, color=MUTED)
  bar(ax, 16, 97, 0.7, TEMPLATE, h=0.2, alpha=0.45)
  for x0 in (85, 77, 69):
    outline(ax, x0, x0 + 8, 0.7, INK, h=0.4, linestyle="--")
  label(ax, 93, 1.12, "the outer probes are not continuous with it",
        fontsize=7, color=OK, ha="right")
  save(fig, "split")


def keep():
  """The part of a read that is kept, in the four cases."""
  fig, ax = canvas(6.6, 3.6, 100, 7.2)
  notes = [(6.3, "Both ends terminal, with references",
            "tails removed,\nconserved termini kept"),
           (4.5, "An end internal or unknown (here the left one), or no "
            "references",
            "at that end the whole\nprimer is removed"),
           (2.7, "Primer dimer", "the primers touch or\noverlap: nothing "
            "is kept"),
           (0.9, "The primer of the same end at both read ends",
            "the second copy, found in\nthe other orientation, is\nremoved "
            "with what lies\nbeyond it")]
  for y, title, note in notes:
    label(ax, 1, y + 0.5, title, fontweight="bold")
    label(ax, 74, y - 0.1, note, fontsize=7, color=MUTED)

  def read(y, left_end=14, right_start=58):
    """Flanks, tails, primers and template of a read, on the row y."""
    bar(ax, 2, left_end - 8, y, FLANK)
    bar(ax, left_end - 8, left_end - 4, y, TAIL)
    bar(ax, left_end - 4, left_end, y, PRIMER)
    bar(ax, right_start, right_start + 4, y, PRIMER)
    bar(ax, right_start + 4, right_start + 8, y, TAIL)
    bar(ax, right_start + 8, 70, y, FLANK)
    if right_start > left_end:
      bar(ax, left_end, right_start, y, TEMPLATE)

  read(6.3)
  bracket(ax, 10, 62, 5.9, "kept", color=OK)
  read(4.5)
  outline(ax, 10, 14, 4.5, BAD, h=0.5)
  bracket(ax, 14, 62, 4.1, "kept", color=OK)
  read(2.7, left_end=36, right_start=36)
  # two amplicons joined without a primer between them, each read from
  # its primer towards the junction
  bar(ax, 2, 6, 0.9, FLANK)
  bar(ax, 6, 10, 0.9, TAIL)
  bar(ax, 10, 14, 0.9, PRIMER)
  bar(ax, 14, 58, 0.9, TEMPLATE)
  ax.plot([36, 36], [0.62, 1.18], color="white", linewidth=1.2)
  bar(ax, 58, 62, 0.9, PRIMER, hatch="////", edgecolor="white")
  bar(ax, 62, 66, 0.9, TAIL)
  bar(ax, 66, 70, 0.9, FLANK)
  label(ax, 60, 1.38, "second copy", fontsize=6.5, color=MUTED,
        ha="center")
  bracket(ax, 10, 58, 0.5, "kept (here with a terminal left end)", color=OK)
  save(fig, "keep")
