# (c) 2026 Giorgio Gonnella
"""Score the read ends classified by flucut against a simulation.

For each simulated read end with a primer, the truth is known
(benchmark/simulate.py): the Uni12 end of an amplicon is terminal (the
start of the segment) or internal (the internal priming site of HA),
the Uni13 end is terminal (the end of the segment). flucut reports the
origin of each end and where it lies on the reference (the tags ol, ml,
or, mr of the trimmed reads; run flucut with --no-split). The true
position is carried over to the reference by a global alignment of the
simulated segment with it (benchmark/pairwise.py).

An end is counted if the simulated read has a primer there and flucut
found it (what happens before the classification does not depend on
the reference). It is correct if its origin is right and its position
is within TOLERANCE bases of the true one, on a reference of the right
segment.

The references may hold several references of a segment (flucut keeps
them all where it cannot choose): an end is then correct if it is
right on the reference it was mapped to. The reference reported for
the segment is the one with the highest identity to the simulated
segment, and the identity is its own.

Usage:
  python3 -I benchmark/score_read_ends.py <simulation dir> \\
      <trimmed reads> <references.fasta>
Output: a TSV line for each segment (segment number, reference, its
identity to the simulated segment, ends, ends with the right origin,
correct ends, ends of unknown origin).
"""

import gzip
import re
import sys

import pairwise

TOLERANCE = 5
UNI12_LEN = 12
UNI13_LEN = 13


def read_fasta(path):
  """(name, sequence, segment number or 0) of each record."""
  records = []
  for line in open(path):
    line = line.strip()
    if line.startswith(">"):
      m = re.search(r"segment (\d)", line)
      records.append([line[1:].split()[0], [], int(m.group(1)) if m else 0])
    elif records:
      records[-1][1].append(line.upper())
  return [(name, "".join(seq), number) for name, seq, number in records]


def true_ends(simdir, genome):
  """For each read, its segment number and the true (side, origin,
  position) of its ends with a primer; the position is that of the
  first base of the segment after the primer (1-based, on the plus
  strand of the simulated segment), as flucut reports it."""
  number = {name: n for name, _, n in genome}
  length = {name: len(seq) for name, seq, _ in genome}
  with open(simdir + "/site.tsv") as f:
    f.readline()
    site = int(f.readline().split("\t")[1])
  ends = {}
  with gzip.open(simdir + "/truth.tsv.gz", "rt") as f:
    f.readline()
    for line in f:
      name, segment, kind, strand, labels = line.rstrip("\n").split("\t")[:5]
      if kind == "dimer" or "T" not in labels:
        continue
      if strand == "-":
        labels = labels[::-1]
      left = labels[:labels.index("T")]
      right = labels[labels.rindex("T"):]
      read_ends = []
      if "i" in left:
        read_ends.append(("l", "internal", site + UNI12_LEN + 1))
      elif "p" in left:
        read_ends.append(("l", "terminal", UNI12_LEN + 1))
      if "p" in right:
        read_ends.append(("r", "terminal", length[segment] - UNI13_LEN))
      ends[name] = (number[segment], read_ends)
  return ends


def read_tags(path):
  """The tags of each trimmed read, by read name."""
  tags = {}
  with gzip.open(path, "rt") as f:
    for n, line in enumerate(f):
      if n % 4 == 0:
        fields = line.rstrip("\n").split("\t")
        tags[fields[0][1:].split()[0]] = {
            t[:2]: t[5:] for t in fields[1:] if len(t) > 5}
  return tags


def score(simdir, trimmed, references, alignments=None):
  """The scores by segment number: a dictionary with reference,
  identity, ends, origin_ok, correct, unknown. <alignments> caches the
  alignments, by (simulated segment, reference)."""
  genome = read_fasta(simdir + "/genome.fasta")
  refs = read_fasta(references)
  if alignments is None:
    alignments = {}
  result, maps = {}, {}
  for name, seq, n in genome:
    candidates = [(r, s) for r, s, rn in refs if rn == n]
    if not candidates:
      sys.exit("no reference of segment %d in %s" % (n, references))
    maps[n] = {}
    best = None
    for ref, ref_seq in candidates:
      if (name, ref) not in alignments:
        alignments[(name, ref)] = pairwise.align(seq, ref_seq)
      identity, maps[n][ref] = alignments[(name, ref)]
      if best is None or identity > best[1]:
        best = (ref, identity)
    result[n] = {"reference": best[0], "identity": best[1], "ends": 0,
                 "origin_ok": 0, "correct": 0, "unknown": 0}
  tags = read_tags(trimmed)
  for name, (n, read_ends) in true_ends(simdir, genome).items():
    t = tags.get(name, {})
    r = result[n]
    for side, origin, position in read_ends:
      if "f" + side not in t:
        continue
      r["ends"] += 1
      found = t.get("o" + side, "unknown")
      if found == "unknown":
        r["unknown"] += 1
      if found != origin:
        continue
      r["origin_ok"] += 1
      mapped = t.get("m" + side, "").rsplit(":", 2)
      if len(mapped) != 3 or mapped[0] not in maps[n]:
        continue
      expected = maps[n][mapped[0]][position - 1] + 1
      if abs(int(mapped[1]) - expected) <= TOLERANCE:
        r["correct"] += 1
  return result


COLUMNS = ["reference", "identity", "ends", "origin_ok", "correct",
           "unknown"]


def main():
  if len(sys.argv) != 4:
    sys.exit(__doc__)
  result = score(*sys.argv[1:])
  print("\t".join(["segment"] + COLUMNS))
  for n in sorted(result):
    r = dict(result[n], identity="%.4f" % result[n]["identity"])
    print("\t".join([str(n)] + [str(r[c]) for c in COLUMNS]))


if __name__ == "__main__":
  main()
