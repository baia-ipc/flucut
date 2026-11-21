# (c) 2026 Giorgio Gonnella
"""Global alignment of two segments: identity and coordinates.

Used by the scripts of the reference choice to compare a simulated
segment with a reference: how similar they are, and which position of
the reference corresponds to a position of the segment (the two differ
by insertions and deletions).

The alignment minimizes the edit distance (substitutions, insertions
and deletions cost 1). Needs numpy.
"""

import numpy as np


def align(a, b):
  """Align the sequences a and b. Returns (identity, positions):
  identity is the share of the alignment columns with equal bases,
  without the gaps at the ends of the alignment (one sequence may lack
  its termini);
  positions[i] is the 0-based position of b aligned to position i of a
  (for a base of a aligned to a gap, the position of b before the
  gap, at least 0)."""
  n, m = len(a), len(b)
  bb = np.frombuffer(b.encode(), dtype=np.uint8)
  index = np.arange(m + 1, dtype=np.int32)
  rows = np.empty((n + 1, m + 1), dtype=np.int32)
  rows[0] = index
  for i in range(1, n + 1):
    prev = rows[i - 1]
    t = np.empty(m + 1, dtype=np.int32)
    t[0] = i
    t[1:] = np.minimum(prev[:-1] + (bb != ord(a[i - 1])), prev[1:] + 1)
    # a gap in a: the best cell to the left, plus the distance to it
    rows[i] = np.minimum.accumulate(t - index) + index
  positions = [0] * n
  # the columns, from the last one: whether the bases are equal, None
  # for a gap
  columns = []
  i, j = n, m
  while i > 0 or j > 0:
    if i > 0 and j > 0 and rows[i, j] == rows[i - 1, j - 1] + (
        a[i - 1] != b[j - 1]):
      columns.append(a[i - 1] == b[j - 1])
      positions[i - 1] = j - 1
      i, j = i - 1, j - 1
    elif i > 0 and rows[i, j] == rows[i - 1, j] + 1:
      columns.append(None)
      positions[i - 1] = max(j - 1, 0)
      i -= 1
    else:
      columns.append(None)
      j -= 1
  # without the gaps at the ends (a sequence without its termini)
  while columns and columns[0] is None:
    columns.pop(0)
  while columns and columns[-1] is None:
    columns.pop()
  return sum(c is True for c in columns) / max(len(columns), 1), positions
