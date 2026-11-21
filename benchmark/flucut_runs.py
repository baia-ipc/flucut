# (c) 2026 Giorgio Gonnella
"""Run the simulator and flucut for the reference-choice benchmarks."""

import gzip
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def simulate(genome, outdir, seed, error_rate, reads=400, segment_reads="",
             prefix=""):
  """Simulate shuffled reads of a genome (benchmark/simulate.py)."""
  command = [sys.executable, "-I", os.path.join(HERE, "simulate.py"),
             genome, outdir, "--seed", str(seed), "--shuffle",
             "--first-site",
             "--error-rate", str(error_rate),
             "--reads-per-segment", str(reads)]
  if segment_reads:
    command += ["--segment-reads", segment_reads]
  if prefix:
    command += ["--read-prefix", prefix]
  subprocess.run(command, check=True)


def mix(outdir, seed, simdirs):
  """Join simulations (benchmark/mix_reads.py)."""
  subprocess.run([sys.executable, "-I", os.path.join(HERE, "mix_reads.py"),
                  outdir, str(seed)] + simdirs, check=True)


def first_reads(reads, n, path):
  """Write the first n reads of a gzipped FASTQ file to path."""
  with gzip.open(reads, "rt") as f, gzip.open(path, "wt") as out:
    for _ in range(4 * n):
      out.write(f.readline())


def flucut(program, reads, references, outdir):
  """Run flucut on the reads against the references; returns its log
  (with the debug messages) and the path of the trimmed reads."""
  result = subprocess.run(
      [program, "-D", "--no-split", "-j", "1", "-r", references,
       "-o", outdir, reads],
      stderr=subprocess.PIPE, stdout=subprocess.DEVNULL, text=True)
  if result.returncode != 0:
    sys.exit("flucut failed on %s:\n%s" % (reads, result.stderr[-2000:]))
  name = re.sub(r"\.fastq\.gz$", "", os.path.basename(reads))
  return result.stderr, os.path.join(outdir, name + ".fastq.gz")


def reference_subtypes(references):
  """The subtype in the description of each reference, by name."""
  subtypes = {}
  for line in open(references):
    if line.startswith(">"):
      m = re.search(r"\((H\d+N\d+)\)", line)
      subtypes[line[1:].split()[0]] = m.group(1) if m else "unknown"
  return subtypes


def split_by_subtype(references, outdir):
  """Write the references of each subtype to a file of their own (one
  reference for each segment); returns the paths by subtype."""
  os.makedirs(outdir, exist_ok=True)
  records = {}
  for record in open(references).read().split(">")[1:]:
    m = re.search(r"\((H\d+N\d+)\)", record.split("\n", 1)[0])
    if m:
      records.setdefault(m.group(1), []).append(">" + record)
  paths = {}
  for subtype, subtype_records in sorted(records.items()):
    paths[subtype] = os.path.join(outdir, subtype + ".fasta")
    with open(paths[subtype], "w") as f:
      f.write("".join(subtype_records))
  return paths
