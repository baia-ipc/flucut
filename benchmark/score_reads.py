# (c) 2026 Giorgio Gonnella
"""Score the trimmed reads of a tool against the truth of a simulation.

Each output read is found by its name (the first word of its title) in
the simulated reads, and located in its input read (an exact substring,
on either strand: the tools cut, they do not change the bases). The
origin of each base of the input read (truth.tsv.gz of
benchmark/simulate.py) then tells what was kept and what was removed:

  reads_in, reads_out   reads simulated, reads written
  reads_lost            reads with template bases which were not written
  template_removed      template bases of the written reads removed
  primer_kept           primer bases kept (tail and terminus, at the
                        segment ends), flank_kept: flanking bases kept
  internal_kept         primer bases kept at the internal priming site
                        (the bases which carry the artefact)
  not_found             written reads not found in their input read

Usage:
  python3 -I benchmark/score_reads.py <simulation dir> <trimmed reads>
      <label> [--no-header]
Output: a TSV line on the standard output (with a header).
"""

import gzip
import sys

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
  return seq.translate(COMPLEMENT)[::-1]


def open_text(path):
  return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fastq(path):
  """(name, sequence) of each record of a FASTQ file"""
  with open_text(path) as f:
    while True:
      title = f.readline()
      if not title:
        return
      seq = f.readline().strip()
      f.readline()
      f.readline()
      yield title[1:].split()[0], seq


COLUMNS = ["reads_in", "reads_out", "reads_lost", "template_removed",
           "primer_kept", "flank_kept", "internal_kept", "not_found"]


def main():
  args = [a for a in sys.argv[1:] if a != "--no-header"]
  if len(args) != 3:
    sys.exit(__doc__)
  simdir, trimmed, label = args
  truth = {}
  with gzip.open(simdir + "/truth.tsv.gz", "rt") as f:
    f.readline()
    for line in f:
      name, _, kind, _, labels = line.rstrip("\n").split("\t")
      truth[name] = (kind, labels)
  inputs = dict(read_fastq(simdir + "/reads.fastq.gz"))
  c = dict.fromkeys(COLUMNS, 0)
  c["reads_in"] = len(truth)
  written = set()
  for name, seq in read_fastq(trimmed):
    c["reads_out"] += 1
    written.add(name)
    read = inputs[name]
    labels = truth[name][1]
    start = read.find(seq)
    if start < 0:
      start = read.find(revcomp(seq))
    if start < 0:
      c["not_found"] += 1
      continue
    end = start + len(seq)
    kept = labels[start:end]
    c["template_removed"] += (labels[:start] + labels[end:]).count("T")
    c["primer_kept"] += kept.count("p") + kept.count("t")
    c["flank_kept"] += kept.count("f")
    c["internal_kept"] += kept.count("i")
  for name, (kind, labels) in truth.items():
    if name not in written and "T" in labels:
      c["reads_lost"] += 1
  if "--no-header" not in sys.argv:
    print("\t".join(["label"] + COLUMNS))
  print("\t".join([label] + [str(c[k]) for k in COLUMNS]))


if __name__ == "__main__":
  main()
