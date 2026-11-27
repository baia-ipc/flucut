# (c) 2026 Giorgio Gonnella
"""Make a synthetic reference of another subtype for the program tests.

The record of <reference.fasta> is copied with a share of its bases
replaced at random (not its first and last KEPT bases, where the
primers bind), under the name and the description given: a sequence
of the same segment which is as far from the original as a reference
of another subtype is. It is not the sequence of a virus; it lets the
tests have reads of two HA subtypes (see testdata/README.md).

The output is deterministic (fixed random seed).

Usage:
  python3 -I scripts/diverged-reference.py <reference.fasta> <name> \
      "<description>"
"""

import random
import sys

SHARE = 0.3
KEPT = 30
WIDTH = 70


def diverged(seq, rng):
  """The sequence with SHARE of its inner bases replaced by another
  base."""
  bases = list(seq)
  for i in range(KEPT, len(bases) - KEPT):
    if rng.random() < SHARE:
      bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
  return "".join(bases)


def main():
  if len(sys.argv) != 4:
    sys.exit(__doc__)
  path, name, description = sys.argv[1:]
  lines = open(path).read().splitlines()
  seq = "".join(l for l in lines[1:] if not l.startswith(">")).upper()
  seq = diverged(seq, random.Random(1))
  print(">%s %s" % (name, description))
  for i in range(0, len(seq), WIDTH):
    print(seq[i:i + WIDTH])


if __name__ == "__main__":
  main()
