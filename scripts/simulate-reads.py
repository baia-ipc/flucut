# (c) 2026 Giorgio Gonnella
"""Simulate the reads of the program tests of flucut.

The reads are built from one segment reference (a real record, see
testdata/README.md) and the primers of the default scheme of flucut
(uni12-mbt13: the untailed Uni12 primer, uni12r in src/primers.c, in
its A form, and the tailed Uni13 primer mbtuni-13). Each read covers
one situation the trimmer has to handle: complete amplicons in both
orientations, flanking bases, sequencing errors, a read starting inside
the primer, internally primed amplicons (by Uni12 and by Uni13, at
arbitrary positions: they are simulated, not the site of a real
artefact), a primer dimer, a read without primers, a short amplicon and
an amplicon with one primer only.

The output is deterministic (fixed random seed); the reads are
simulated, not sequenced.

Usage:
  python3 -I scripts/simulate-reads.py <reference.fasta> > reads.fastq
"""

import random
import sys

# the primers of the default scheme, as found in a read of the plus
# strand (src/primers.h): Uni12 at the start, the reverse complement of
# the tailed Uni13 primer at the end
UNI12 = "AGCAAAAGCAGG"
RC_MBT13 = "CCTTGTTTCTACTGATCACGCGT"
# the conserved terminus of the Uni12 end
TERMINUS12 = "AGCAAAAGCAGG"

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
  """Reverse complement of a DNA sequence."""
  return seq.translate(COMPLEMENT)[::-1]


def read_fasta(path):
  """The sequence of the first record of a FASTA file."""
  seq = []
  with open(path) as f:
    for i, line in enumerate(f):
      if line.startswith(">"):
        if i > 0:
          break
        continue
      seq.append(line.strip().upper())
  return "".join(seq)


def template(ref):
  """The reference without the part of the Uni12 terminus it contains
  at its start (the primer supplies it in a read)."""
  for k in range(len(TERMINUS12), 0, -1):
    if ref.startswith(TERMINUS12[-k:]):
      return ref[k:]
  return ref


def mutate(seq, rate, rng):
  """Substitutions at the given rate."""
  out = list(seq)
  for i, c in enumerate(out):
    if rng.random() < rate:
      out[i] = rng.choice([b for b in "ACGT" if b != c])
  return "".join(out)


def main():
  if len(sys.argv) != 2:
    sys.exit(__doc__)
  rng = random.Random(12)
  t = template(read_fasta(sys.argv[1]))

  def flank(n):
    return "".join(rng.choice("ACGT") for _ in range(n))

  full = UNI12 + t + RC_MBT13
  internal12 = UNI12 + t[700:] + RC_MBT13
  reads = [
      ("full_plus", full),
      ("full_minus", revcomp(full)),
      ("flanked_plus", flank(15) + full + flank(10)),
      ("errors_plus", mutate(full, 0.04, rng)),
      ("partial_left", UNI12[5:] + t + RC_MBT13),
      ("internal12_plus", internal12),
      ("internal12_minus", flank(12) + revcomp(internal12)),
      ("internal13_plus", UNI12 + t[:900] + RC_MBT13),
      ("dimer", UNI12 + RC_MBT13),
      ("no_primer", t[100:400]),
      ("short_amplicon", UNI12 + t[:20] + RC_MBT13),
      ("left_only", UNI12 + t[:500]),
  ]
  for name, seq in reads:
    qual = "".join(chr(33 + 5 + (i * 7) % 35) for i in range(len(seq)))
    sys.stdout.write("@%s simulated\n%s\n+\n%s\n" % (name, seq, qual))


if __name__ == "__main__":
  main()
