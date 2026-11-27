# (c) 2026 Giorgio Gonnella
"""Measure what flucut reports about the reads of a second strain in
mixed samples.

flucut reports, for the HA and the NA, the reads of a second subtype
among the first reads of a sample (src/subtype-reads.h; HTML report and
JSON summary), and, for each segment, how often the probes of one read
vote for the same reference (src/read-agreement.h; messages only). Here
samples of one strain and mixes of two strains are simulated, to see
how well the two statistics tell them apart.

Each sample has <reads> reads per segment in all (default 400); in a
mix, a share of them is of the minor strain (MINOR_FRACTIONS). With
<NA reads>, the NA has that number of reads instead (a sample whose NA
is scarce: flucut then counts the probes of up to 20,000 reads instead
of 2000). All genomes of <genomes dir> are simulated alone; the mixes
are those of PAIRS.

Output, in <workdir>:
- second_subtype.tsv, with a line for each sample and for the HA and
  the NA, from the JSON summary: major and minor strain (minor: none
  for one strain), minor fraction, error rate, seed, gene, the subtype
  of the minor strain (none: no second subtype in the sample), then
  what flucut reports: leading and second subtype, reads with several
  voting probes, reads of the leading and of the second subtype, those
  expected, the excess, its share of the reads, and whether it is
  above the threshold;
- agreement.tsv, with a line for each sample, segment and reference,
  from the messages: the same sample columns, segment, reads with
  several voting probes, reference, its subtype, its votes, the votes
  of the segment, the reads agreeing on it, and the reads expected to
  if the probes voted independently of their read.

Usage:
  python3 benchmark/mixed_samples.py <genomes dir> <references.fasta> \\
      <flucut> <workdir> [<jobs> [<reads> [<NA reads>]]]
"""

import glob
import json
import multiprocessing
import os
import re
import shutil
import sys
import zlib

import flucut_log
import flucut_runs

ERROR_RATES = [0.04, 0.08, 0.12]
SEEDS = [1, 2]
MINOR_FRACTIONS = [0.5, 0.1, 0.05, 0.02, 0.01]
GENES = {"ha": "H", "na": "N"}
# (major, minor): strains of different subtypes, except the last pair
# (two lineages of H5N1); the swine pair shares its internal genes
PAIRS = [("ny04_h3n2", "vn12_h5n1"), ("vn12_h5n1", "ny04_h3n2"),
         ("vn12_h5n1", "km17_h7n9"), ("km17_h7n9", "vn13_h5n1"),
         ("ah15_h9n2", "vn12_h5n1"), ("km17_h7n9", "ah15_h9n2"),
         ("pr8_h1n1", "ny04_h3n2"), ("cl14_h1n2", "mx12_h3n2"),
         ("vn12_h5n1", "vn13_h5n1")]


def true_second(major, minor, letter):
  """The subtype of the minor strain in the HA (letter H) or the NA (N),
  from the labels of the genomes (e.g. vn12_h5n1); none if the sample
  has one strain or both strains have the same."""
  def number(label):
    return re.search(letter + r"(\d+)", label.split("_")[1].upper()).group(0)
  if minor == "none" or number(major) == number(minor):
    return "none"
  return number(minor)


def subtype_rows(sample, outdir):
  """The rows of second_subtype.tsv of a sample, from its JSON summary."""
  major, minor = sample[0], sample[1]
  with open(os.path.join(outdir, "flucut_summary.json")) as f:
    reads = json.load(f)["all"]["sample"]["subtype_reads"]
  rows = []
  for gene, letter in GENES.items():
    e = reads[gene]
    truth = true_second(major, minor, letter)
    if e is None:
      rows.append(sample + [gene.upper(), truth] + ["NA"] * 9)
      continue
    rows.append(sample + [gene.upper(), truth, e["leading"],
                          e["second"] or "none", e["reads"],
                          e["leading_reads"], e["second_reads"],
                          e["second_expected"], e["second_excess"],
                          e["second_share"], int(e["above_threshold"])])
  return rows


def segment_reads(n_reads, n_na, total_reads):
  """The --segment-reads value which gives the NA its share of n_na
  reads in total_reads per segment ("": the NA as the other segments)."""
  if n_na is None:
    return ""
  return "6:%d" % round(n_na * n_reads / total_reads)


def run(task):
  """The rows of second_subtype.tsv and of agreement.tsv of a sample."""
  (major, minor, fraction, error_rate, seed, genomes, references, program,
   workdir, n_reads, n_na) = task
  subtypes = flucut_runs.reference_subtypes(references)
  base = os.path.join(workdir, "samples", "%s_%s_%g_%g_%d" % (
      major, minor, fraction, error_rate, seed))
  n_minor = round(n_reads * fraction) if minor != "none" else 0
  flucut_runs.simulate("%s/%s.fasta" % (genomes, major), base + "/major",
                       seed, error_rate, reads=n_reads - n_minor,
                       segment_reads=segment_reads(n_reads - n_minor, n_na,
                                                   n_reads), prefix="a")
  sample = base + "/major"
  if n_minor > 0:
    flucut_runs.simulate("%s/%s.fasta" % (genomes, minor), base + "/minor",
                         seed + 1000, error_rate, reads=n_minor,
                         segment_reads=segment_reads(n_minor, n_na, n_reads),
                         prefix="b")
    sample = base + "/mix"
    # (an order of its own for each sample: with the seed alone, the
    # reads of the minor strain would be at the same places in every
    # mix of one seed, and flucut counts the first reads only)
    flucut_runs.mix(sample, zlib.crc32(os.path.basename(base).encode()),
                    [base + "/major", base + "/minor"])
  log, _ = flucut_runs.flucut(program, sample + "/reads.fastq.gz",
                              references, base + "/out")
  sample = [major, minor, fraction if n_minor else 0, error_rate, seed]
  rows = []
  for segment, e in sorted(flucut_log.parse(log).items()):
    for ref, votes in e["votes"].items():
      rows.append(sample + [segment, e["reads"], ref, subtypes[ref], votes,
                            sum(e["votes"].values()),
                            e["agreeing"].get(ref, 0),
                            "%.1f" % e["expected"].get(ref, 0)])
  by_subtype = subtype_rows(sample, base + "/out")
  shutil.rmtree(base)
  return by_subtype, rows


SAMPLE_COLUMNS = ["major", "minor", "minor_fraction", "error_rate", "seed"]


def write_table(path, header, rows):
  with open(path, "w") as f:
    f.write("\t".join(header) + "\n")
    for row in rows:
      f.write("\t".join(str(v) for v in row) + "\n")


def main():
  if not 5 <= len(sys.argv) <= 8:
    sys.exit(__doc__)
  genomes, references, program, workdir = sys.argv[1:5]
  jobs = int(sys.argv[5]) if len(sys.argv) > 5 else 4
  n_reads = int(sys.argv[6]) if len(sys.argv) > 6 else 400
  n_na = int(sys.argv[7]) if len(sys.argv) > 7 else None
  labels = [os.path.basename(g)[:-len(".fasta")]
            for g in sorted(glob.glob(genomes + "/*.fasta"))]
  samples = [(label, "none", 0) for label in labels]
  samples += [(major, minor, f) for major, minor in PAIRS
              for f in MINOR_FRACTIONS]
  tasks = [(major, minor, f, e, seed, genomes, references, program, workdir,
            n_reads, n_na)
           for major, minor, f in samples for e in ERROR_RATES
           for seed in SEEDS]
  os.makedirs(workdir, exist_ok=True)
  with multiprocessing.Pool(jobs) as pool:
    results = pool.map(run, tasks, chunksize=1)
  write_table(os.path.join(workdir, "second_subtype.tsv"),
              SAMPLE_COLUMNS + ["gene", "minor_subtype", "leading", "second",
                                "reads", "leading_reads", "second_reads",
                                "second_expected", "second_excess",
                                "second_share", "above_threshold"],
              [row for by_subtype, _ in results for row in by_subtype])
  write_table(os.path.join(workdir, "agreement.tsv"),
              SAMPLE_COLUMNS + ["segment", "reads", "reference", "subtype",
                                "votes", "segment_votes", "agreeing",
                                "expected"],
              [row for _, rows in results for row in rows])


if __name__ == "__main__":
  main()
