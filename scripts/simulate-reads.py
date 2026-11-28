# (c) 2026 Giorgio Gonnella
"""Simulate the reads of the program tests of flucut.

The reads are built from one segment reference (a real record, see
testdata/README.md) and one pair of the primers built into flucut
(src/seq/primers/primers.c). Each read covers one situation the trimmer has to
handle: complete amplicons in both orientations, flanking bases,
sequencing errors, a read starting inside the primer, internally primed
amplicons (at arbitrary positions: they are simulated, not the site of
a real artefact), a primer dimer, a read without primers, a short
amplicon and an amplicon with one primer only.

The default pair is that of the default scheme of influenza A
(uni12-mbt13: the untailed Uni12 primer, uni12r in its A form, and the
tailed Uni13 primer mbtuni-13). --primers chooses a pair of influenza B
(the IBV-GA2 cocktail) instead; the pairs of influenza B each bind
their own segments, so the pair must be that of the reference. For them
one further read is simulated, primed by the pair of another segment:
its primer bases are not the sequence of this segment, and flucut has
to remove them although the read starts at a segment terminus.

With the argument "deletion", the reads are instead amplicons with an
internal deletion (the two ends of the segment joined, as in a
defective interfering genome), in both orientations.

The output is deterministic (fixed random seed); the reads are
simulated, not sequenced.

Usage:
  python3 -I scripts/simulate-reads.py [--primers <pair>] \
      <reference.fasta> [deletion]

  <pair>: uni12-mbt13 (default), b-pbs, b-pa, b-hana, b-np, b-m or b-ns
"""

import random
import sys

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
  """Reverse complement of a DNA sequence."""
  return seq.translate(COMPLEMENT)[::-1]


# The pairs of primers, as they are ordered (5' -> 3'), with the number
# of their trailing bases which are the sequence of the segment (the
# rest is the tail of the protocol); they are those built into flucut,
# see src/seq/primers/primers.c. A read of the plus strand of an amplicon begins
# with the left primer and ends with the reverse complement of the
# right one.
PAIRS = {
  "uni12-mbt13": (("AGCAAAAGCAGG", 12), ("ACGCGTGATCAGTAGAAACAAGG", 13)),
  "b-pbs": (("GGGGGGAGCAGAAGCGGAGC", 14),
            ("CCGGGTTATTAGTAGAAACACGAGC", 15)),
  "b-pa": (("GGGGGGAGCAGAAGCGGTGC", 14),
           ("CCGGGTTATTAGTAGAAACACGTGC", 15)),
  "b-hana": (("GGGGGGAGCAGAAGCAGAGC", 14),
             ("CCGGGTTATTAGTAGTAACAAGAGC", 15)),
  "b-np": (("GGGGGGAGCAGAAGCACAGC", 14),
           ("CCGGGTTATTAGTAGAAACAACAGC", 15)),
  "b-m": (("GGGGGGAGCAGAAGCACGCACTT", 17),
          ("CCGGGTTATTAGTAGAAACAACGCACTT", 18)),
  "b-ns": (("GGGGGGAGCAGAAGCAGAGGATT", 17),
           ("CCGGGTTATTAGTAGTAACAAGAGGATT", 18)),
}

# the pair used for the read primed by the primers of another segment
OTHER_PAIR = {"b-hana": "b-m", "b-m": "b-hana", "b-pbs": "b-ns",
              "b-pa": "b-np", "b-np": "b-pa", "b-ns": "b-pbs"}


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


# the fewest bases of a terminus recognised at an end of a reference:
# fewer could be there by chance, and would cut a base of the template
MIN_TERMINUS = 4


def template(ref, left, right):
  """The reference without the bases of the segment termini it carries,
  at its start and at its end: the primers supply them in a read. A
  record may carry all of them, some of them or none."""
  for k in range(len(left), MIN_TERMINUS - 1, -1):
    if ref.startswith(left[-k:]):
      ref = ref[k:]
      break
  for k in range(len(right), MIN_TERMINUS - 1, -1):
    if ref.endswith(right[:k]):
      ref = ref[:-k]
      break
  return ref


def mutate(seq, rate, rng):
  """Substitutions at the given rate."""
  out = list(seq)
  for i, c in enumerate(out):
    if rng.random() < rate:
      out[i] = rng.choice([b for b in "ACGT" if b != c])
  return "".join(out)


def deletion_reads(t, left, right):
  """Amplicons of the template <t> without its middle part."""
  deleted = left + t[:250] + t[-300:] + right
  return [
      ("deletion_plus", deleted),
      ("deletion_minus", revcomp(deleted)),
  ]


def situation_reads(t, left, right, other_left):
  """A read for each situation of the trimming of the template <t>.

  <left> and <right> are the primers as they are found in a read of the
  plus strand; <other_left>, where there is one, the left primer of the
  segments of another pair."""
  rng = random.Random(12)

  def flank(n):
    return "".join(rng.choice("ACGT") for _ in range(n))

  full = left + t + right
  internal_left = left + t[700:] + right
  reads = [
      ("full_plus", full),
      ("full_minus", revcomp(full)),
      ("flanked_plus", flank(15) + full + flank(10)),
      ("errors_plus", mutate(full, 0.04, rng)),
      ("partial_left", left[5:] + t + right),
      ("internal12_plus", internal_left),
      ("internal12_minus", flank(12) + revcomp(internal_left)),
      ("internal13_plus", left + t[:900] + right),
      ("dimer", left + right),
      ("no_primer", t[100:400]),
      ("short_amplicon", left + t[:20] + right),
      ("left_only", left + t[:500]),
  ]
  if other_left is not None:
    # the amplicon began at this segment terminus, but was primed by
    # the primer of other segments: its bases are not this segment
    reads.append(("other_segment_primer", other_left + t + right))
  return reads


def parse_args(argv):
  """The pair of primers, the reference and whether the reads are the
  amplicons with an internal deletion."""
  pair = "uni12-mbt13"
  if argv[:1] == ["--primers"]:
    if len(argv) < 2 or argv[1] not in PAIRS:
      sys.exit(__doc__)
    pair, argv = argv[1], argv[2:]
  if len(argv) not in (1, 2) or argv[1:] not in ([], ["deletion"]):
    sys.exit(__doc__)
  return pair, argv[0], len(argv) == 2


def main():
  pair, path, deletion = parse_args(sys.argv[1:])
  (left, left_len), (right, right_len) = PAIRS[pair]
  right = revcomp(right)
  # the bases of the segment in each primer, as they are in the plus
  # strand: the suffix of the left primer, the prefix of the right one
  # once it is reverse-complemented
  t = template(read_fasta(path), left[-left_len:], right[:right_len])
  other = PAIRS[OTHER_PAIR[pair]][0][0] if pair in OTHER_PAIR else None
  reads = deletion_reads(t, left, right) if deletion else \
    situation_reads(t, left, right, other)
  for name, seq in reads:
    qual = "".join(chr(33 + 5 + (i * 7) % 35) for i in range(len(seq)))
    sys.stdout.write("@%s simulated\n%s\n+\n%s\n" % (name, seq, qual))


if __name__ == "__main__":
  main()
