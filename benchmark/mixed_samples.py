# (c) 2026 Giorgio Gonnella
"""Measure the agreement of the probes within a read in mixed samples.

flucut reports, for each segment, how often the probes of one read vote
for the same reference (src/read-agreement.h). Here samples of one
strain and mixes of two strains are simulated, to see whether that
statistic tells them apart: a split vote of one strain lying between
two references, against the split vote of two strains.

Each sample has 400 reads per segment in all; in a mix, a share of
them is of the minor strain (MINOR_FRACTIONS). All genomes of
<genomes dir> are simulated alone; the mixes are those of PAIRS.

Output: <workdir>/agreement.tsv, with a line for each sample, segment
and reference: major and minor strain (minor: none for one strain),
minor fraction, error rate, seed, segment, reads with several voting
probes, reference, its subtype, its votes, the votes of the segment,
the reads agreeing on it, and the reads expected to if the probes voted
independently of their read.

Usage:
  python3 benchmark/mixed_samples.py <genomes dir> <references.fasta> \\
      <flucut> <workdir> [<jobs>]
"""

import glob
import multiprocessing
import os
import shutil
import sys

import flucut_log
import flucut_runs

ERROR_RATES = [0.04, 0.08, 0.12]
SEEDS = [1, 2]
MINOR_FRACTIONS = [0.5, 0.1, 0.05, 0.02, 0.01]
READS = 400
# (major, minor): strains of different subtypes, except the last pair
# (two lineages of H5N1); the swine pair shares its internal genes
PAIRS = [("ny04_h3n2", "vn12_h5n1"), ("vn12_h5n1", "ny04_h3n2"),
         ("vn12_h5n1", "km17_h7n9"), ("km17_h7n9", "vn13_h5n1"),
         ("ah15_h9n2", "vn12_h5n1"), ("km17_h7n9", "ah15_h9n2"),
         ("pr8_h1n1", "ny04_h3n2"), ("cl14_h1n2", "mx12_h3n2"),
         ("vn12_h5n1", "vn13_h5n1")]


def run(task):
  """The rows of agreement.tsv of a sample."""
  (major, minor, fraction, error_rate, seed, genomes, references, program,
   workdir) = task
  subtypes = flucut_runs.reference_subtypes(references)
  base = os.path.join(workdir, "samples", "%s_%s_%g_%g_%d" % (
      major, minor, fraction, error_rate, seed))
  n_minor = round(READS * fraction) if minor != "none" else 0
  flucut_runs.simulate("%s/%s.fasta" % (genomes, major), base + "/major",
                       seed, error_rate, reads=READS - n_minor, prefix="a")
  sample = base + "/major"
  if n_minor > 0:
    flucut_runs.simulate("%s/%s.fasta" % (genomes, minor), base + "/minor",
                         seed + 1000, error_rate, reads=n_minor, prefix="b")
    sample = base + "/mix"
    flucut_runs.mix(sample, seed, [base + "/major", base + "/minor"])
  log, _ = flucut_runs.flucut(program, sample + "/reads.fastq.gz",
                              references, base + "/out")
  rows = []
  for segment, e in sorted(flucut_log.parse(log).items()):
    for ref, votes in e["votes"].items():
      rows.append([major, minor, fraction if n_minor else 0, error_rate,
                   seed, segment, e["reads"], ref, subtypes[ref], votes,
                   sum(e["votes"].values()), e["agreeing"].get(ref, 0),
                   "%.1f" % e["expected"].get(ref, 0)])
  shutil.rmtree(base)
  return rows


def main():
  if len(sys.argv) not in (5, 6):
    sys.exit(__doc__)
  genomes, references, program, workdir = sys.argv[1:5]
  jobs = int(sys.argv[5]) if len(sys.argv) == 6 else 4
  labels = [os.path.basename(g)[:-len(".fasta")]
            for g in sorted(glob.glob(genomes + "/*.fasta"))]
  samples = [(label, "none", 0) for label in labels]
  samples += [(major, minor, f) for major, minor in PAIRS
              for f in MINOR_FRACTIONS]
  tasks = [(major, minor, f, e, seed, genomes, references, program, workdir)
           for major, minor, f in samples for e in ERROR_RATES
           for seed in SEEDS]
  os.makedirs(workdir, exist_ok=True)
  with multiprocessing.Pool(jobs) as pool:
    results = pool.map(run, tasks, chunksize=1)
  with open(os.path.join(workdir, "agreement.tsv"), "w") as f:
    f.write("\t".join(["major", "minor", "minor_fraction", "error_rate",
                       "seed", "segment", "reads", "reference", "subtype",
                       "votes", "segment_votes", "agreeing", "expected"])
            + "\n")
    for rows in results:
      for row in rows:
        f.write("\t".join(str(v) for v in row) + "\n")


if __name__ == "__main__":
  main()
