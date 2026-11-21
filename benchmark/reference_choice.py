# (c) 2026 Giorgio Gonnella
"""Measure what the choice of the reference of each segment costs.

flucut chooses a reference for each segment from the votes of probes of
the first reads (src/auto-settings.h). Here reads are simulated from
genomes which are not in the reference set, so that the truth is known,
and two things are measured for each genome and error rate:

- the accuracy with each reference: the reads (400 for each segment)
  are classified against the references of one subtype at a time, and
  the read ends are scored (benchmark/score_read_ends.py); the whole
  set is scored as well (the subtype "all": what a segment whose
  reference is not chosen gets, as flucut then keeps all of them);
- the votes: flucut is run against the whole set, on samples with
  different numbers of reads (all segments alike; few reads of some
  segments; the first reads only), and the votes of every reference and
  the reference chosen are taken from its log.

Outputs, in <workdir>:
  accuracy.tsv   genome, error rate, segment, subtype, reference, its
                 identity to the simulated segment, ends, ends with the
                 right origin, correct ends, ends of unknown origin
  votes.tsv      genome, error rate, sample, segment, reference, its
                 subtype, its votes, the votes of the segment, the
                 probes which match several references, whether it was
                 chosen, and how
The simulated reads are written below <workdir>/samples and removed
when scored.

Usage:
  python3 benchmark/reference_choice.py <genomes dir> <references.fasta> \\
      <flucut> <workdir> [<jobs>]
"""

import glob
import multiprocessing
import os
import shutil
import sys

import flucut_log
import flucut_runs
import score_read_ends

ERROR_RATES = [0.04, 0.08, 0.12]
# samples with few reads of some segments (see simulate.py)
UNEVEN = {"few_polymerase": "1:6,2:12,3:25,5:50",
          "few_surface": "4:6,6:12,7:25,8:50"}
FIRST_READS = [30, 60, 125, 250, 500, 1000]


def vote_rows(key, sample, log, subtypes):
  rows = []
  for segment, e in sorted(flucut_log.parse(log).items()):
    for ref, votes in e["votes"].items():
      rows.append(key + [sample, segment, ref, subtypes[ref], votes,
                         e["total"], e["ambiguous"],
                         int(ref == e["chosen"]), e["how"]])
  return rows


def run(task):
  """The rows of accuracy.tsv and votes.tsv of a genome and error
  rate."""
  genome, error_rate, references, by_subtype, program, workdir = task
  label = os.path.basename(genome)[:-len(".fasta")]
  key = [label, error_rate]
  subtypes = flucut_runs.reference_subtypes(references)
  base = os.path.join(workdir, "samples", "%s_%g" % (label, error_rate))
  seed = 1 + int(100 * error_rate)
  accuracy, votes, alignments = [], [], {}
  sim = base + "/even"
  flucut_runs.simulate(genome, sim, seed, error_rate)
  reads = sim + "/reads.fastq.gz"
  for subtype, path in list(by_subtype.items()) + [("all", references)]:
    _, trimmed = flucut_runs.flucut(program, reads, path, sim + "/" + subtype)
    scores = score_read_ends.score(sim, trimmed, path, alignments)
    for segment, s in sorted(scores.items()):
      accuracy.append(key + [segment, subtype, s["reference"],
                             "%.4f" % s["identity"], s["ends"],
                             s["origin_ok"], s["correct"], s["unknown"]])
  log, _ = flucut_runs.flucut(program, reads, references, sim + "/all")
  votes += vote_rows(key, "even", log, subtypes)
  for n in FIRST_READS:
    part = "%s/first_%d.fastq.gz" % (sim, n)
    flucut_runs.first_reads(reads, n, part)
    log, _ = flucut_runs.flucut(program, part, references, sim + "/part")
    votes += vote_rows(key, "first_%d" % n, log, subtypes)
  for sample, segment_reads in UNEVEN.items():
    sim = base + "/" + sample
    flucut_runs.simulate(genome, sim, seed, error_rate,
                         segment_reads=segment_reads)
    log, _ = flucut_runs.flucut(program, sim + "/reads.fastq.gz",
                                references, sim + "/all")
    votes += vote_rows(key, sample, log, subtypes)
  shutil.rmtree(base)
  return accuracy, votes


def write(path, header, rows):
  with open(path, "w") as f:
    f.write("\t".join(header) + "\n")
    for row in rows:
      f.write("\t".join(str(v) for v in row) + "\n")


def main():
  if len(sys.argv) not in (5, 6):
    sys.exit(__doc__)
  genomes, references, program, workdir = sys.argv[1:5]
  jobs = int(sys.argv[5]) if len(sys.argv) == 6 else 4
  by_subtype = flucut_runs.split_by_subtype(
      references, os.path.join(workdir, "references"))
  tasks = [(g, e, references, by_subtype, program, workdir)
           for g in sorted(glob.glob(genomes + "/*.fasta"))
           for e in ERROR_RATES]
  with multiprocessing.Pool(jobs) as pool:
    results = pool.map(run, tasks, chunksize=1)
  write(os.path.join(workdir, "accuracy.tsv"),
        ["genome", "error_rate", "segment", "subtype", "reference",
         "identity", "ends", "origin_ok", "correct", "unknown"],
        [row for accuracy, _ in results for row in accuracy])
  write(os.path.join(workdir, "votes.tsv"),
        ["genome", "error_rate", "sample", "segment", "reference",
         "subtype", "votes", "segment_votes", "ambiguous", "chosen", "how"],
        [row for _, votes in results for row in votes])


if __name__ == "__main__":
  main()
