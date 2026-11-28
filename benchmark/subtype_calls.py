# (c) 2026 Giorgio Gonnella
"""Measure whether flucut names the subtype of a sample correctly.

flucut names the subtype of a sample after the references it chooses
for the HA and the NA (src/core/auto-settings.h). A sample of a subtype,
lineage or region which the built-in references do not cover is named
after the closest reference in them, which can be another subtype. This
script measures how often that happens: reads are simulated from the
genomes of a table (benchmark/reference_choice_genomes.tsv: subtypes,
hosts, regions and decades, none of them in the reference set), flucut
is run against the references, and the subtype it reports is compared
with the subtype of the genome.

The HA and the NA are counted separately: a genome is named correctly
only if both are right.

Two numbers of the run tell whether the references also place the read
ends, which is what they are there for: the share of the ends whose
origin flucut cannot determine, and the share of the reads of the
segments other than the HA which it takes for internally primed. The
simulation puts an internal priming site in the HA only, so every
internal end reported elsewhere is one the references got wrong.

Outputs, in <workdir>:
  subtype_calls.tsv  genome, strain, its subtype, the error rate, the
                     subtype reported, the references chosen for the HA
                     and the NA and their subtypes, whether the H and
                     the N are right, the ends not located and the
                     internal priming reported outside the HA
The simulated reads are written below <workdir>/samples and removed
once the sample is named.

Usage:
  python3 benchmark/subtype_calls.py <genomes dir> <genomes.tsv> \\
      <references.fasta> <flucut> <workdir> [<jobs>] [<reads>]
"""

import json
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import zlib

import flucut_runs

ERROR_RATES = [0.04, 0.08, 0.12]
READS_PER_SEGMENT = 400
SUBTYPE = re.compile(r"(H\d+)(N\d+)")


def genomes_of(path):
  """The label, strain and subtype of each genome of the table."""
  out = []
  for line in open(path):
    if line.startswith("#") or line.startswith("label\t"):
      continue
    fields = line.rstrip("\n").split("\t")
    if len(fields) >= 3:
      out.append((fields[0], fields[1], fields[2]))
  return out


def reference_subtypes(path):
  """The subtype in the description of each reference of the FASTA."""
  out = {}
  for line in open(path):
    if line.startswith(">"):
      name = line[1:].split()[0]
      found = SUBTYPE.search(line)
      out[name] = found.group(0) if found else ""
  return out


def run_flucut(program, reads, references, outdir):
  """flucut on the reads against the references; its sample summary."""
  done = subprocess.run([program, "--no-split", "-j", "1", "-r", references,
                         "-o", outdir, reads], stdout=subprocess.DEVNULL,
                        stderr=subprocess.PIPE, text=True)
  if done.returncode != 0:
    sys.exit("flucut failed on %s:\n%s" % (reads, done.stderr[-2000:]))
  with open(os.path.join(outdir, "flucut_summary.json")) as handle:
    return json.load(handle)["all"]["sample"]


def segment_reference(summary, gene):
  for segment in summary.get("segments", []):
    if segment.get("gene", "").startswith(gene):
      return segment.get("reference", "")
  return ""


def false_internal(summary):
  """The share of the reads of the segments other than the HA which
  flucut reports as internally primed: the simulation primes inside the
  HA only, so these are read ends the references place wrongly. A
  segment flucut did not find has no "internal" share and is left out,
  rather than counted as 0."""
  shares = [s["internal"] for s in summary.get("segments", [])
            if not s.get("gene", "").startswith("HA") and "internal" in s]
  return round(sum(shares) / len(shares), 4) if shares else ""


def call_genome(task):
  """Simulate one sample and name it; one row of the output."""
  (label, strain, subtype, error_rate, genomes, references, program,
   workdir, reads_per_segment) = task
  folder = os.path.join(workdir, "samples", "%s_%g" % (label, error_rate))
  os.makedirs(folder, exist_ok=True)
  flucut_runs.simulate(os.path.join(genomes, label + ".fasta"), folder,
                       seed=zlib.crc32(label.encode()) % 10000,
                       error_rate=error_rate, reads=reads_per_segment)
  summary = run_flucut(program, os.path.join(folder, "reads.fastq.gz"),
                       references, os.path.join(folder, "flucut"))
  shutil.rmtree(folder, ignore_errors=True)
  subtypes = reference_subtypes(references)
  ha = segment_reference(summary, "HA")
  na = segment_reference(summary, "NA")
  called = summary.get("subtype", "")
  truth = SUBTYPE.search(subtype)
  h_ok = n_ok = ""
  if truth and called:
    found = SUBTYPE.search(called)
    h_ok = int(bool(found) and found.group(1) == truth.group(1))
    n_ok = int(bool(found) and found.group(2) == truth.group(2))
  return [label, strain, subtype, error_rate, called, ha,
          subtypes.get(ha, ""), na, subtypes.get(na, ""), h_ok, n_ok,
          summary.get("ends_not_located", ""), false_internal(summary)]


def main():
  if len(sys.argv) < 6:
    sys.exit(__doc__)
  genomes, table, references, program, workdir = sys.argv[1:6]
  jobs = int(sys.argv[6]) if len(sys.argv) > 6 else 4
  reads = int(sys.argv[7]) if len(sys.argv) > 7 else READS_PER_SEGMENT
  os.makedirs(workdir, exist_ok=True)
  tasks = [(label, strain, subtype, rate, genomes, references, program,
            workdir, reads)
           for label, strain, subtype in genomes_of(table)
           for rate in ERROR_RATES]
  with multiprocessing.Pool(jobs) as pool:
    rows = pool.map(call_genome, tasks)
  out = os.path.join(workdir, "subtype_calls.tsv")
  with open(out, "w") as handle:
    handle.write("genome\tstrain\tsubtype\terror_rate\tcalled\tha_ref\t"
                 "ha_ref_subtype\tna_ref\tna_ref_subtype\th_correct\t"
                 "n_correct\tends_not_located\tfalse_internal\n")
    for row in rows:
      handle.write("\t".join(str(v) for v in row) + "\n")
  shutil.rmtree(os.path.join(workdir, "samples"), ignore_errors=True)
  both = sum(1 for r in rows if r[9] == 1 and r[10] == 1)
  located = [r[11] for r in rows if r[11] != ""]
  internal = [r[12] for r in rows if r[12] != ""]
  print("%d of %d samples named correctly (H and N); ends not located "
        "%.3f, internal priming outside the HA %.3f (means); %s" %
        (both, len(rows), sum(located) / len(located) if located else 0,
         sum(internal) / len(internal) if internal else 0, out))


if __name__ == "__main__":
  main()
