# (c) 2026 Giorgio Gonnella
"""Compare the consensus of an IRMA run with the truth of a simulation.

For each segment, the amended consensus of IRMA
(amended_consensus/<sample>_<n>.fa, n = segment number), or with
--plain the majority consensus (A_<gene>*.fasta), is aligned to the
true segment (genome.fasta of benchmark/simulate.py, whose headers
contain "segment <n>"; N in the consensus matches any base), and the
following are reported:

  missing_5, missing_3   true bases missing at each end of the consensus
  extra_5, extra_3       consensus bases beyond the ends of the truth
  subst, ins, del        differences inside the aligned part
  ambiguous              consensus bases which are not A, C, G, T or N
                         (counted apart, not as substitutions)
  n_bases                N in the consensus (e.g. the amended consensus
                         where fewer reads than MIN_CONS_SUPPORT of IRMA
                         support a base: 50 with FLU-minion)
  site_errors            for the segment of the internal priming site
                         (site.tsv): differences in the replaced bases
                         (NA for the other segments)

Usage:
  python3 -I benchmark/compare_consensus.py [--plain] <simulation dir> \\
      <IRMA output dir> <label>
Output: a TSV table on the standard output (one line per segment, with
the label in the first column; a header unless --no-header).
"""

import glob
import os
import re
import sys

import edlib


def read_fasta(path):
  """The records of a FASTA file, as (header, sequence) pairs."""
  records = []
  for line in open(path):
    line = line.strip()
    if line.startswith(">"):
      records.append([line[1:], []])
    elif records:
      records[-1][1].append(line.upper())
  return [(h, "".join(s)) for h, s in records]


def segment_number(header):
  m = re.search(r"segment (\d)", header)
  return int(m.group(1)) if m else None


def compare(truth, cons, site):
  """The measures of one segment (see the module documentation)."""
  # the consensus, end-free in the truth, then the ends beyond it
  eq = [("N", b) for b in "ACGT"]
  res = edlib.align(cons, truth, mode="HW", task="path",
                    additionalEqualities=eq)
  start, end = res["locations"][0]
  aln = edlib.getNiceAlignment(res, cons, truth)
  q, t = aln["query_aligned"], aln["target_aligned"]
  # (with HW, the consensus is aligned whole; bases of the consensus
  # aligned to gaps at the ends of the alignment are beyond the truth)
  lead = 0
  while lead < len(t) and t[lead] == "-":
    lead += 1
  trail = 0
  while trail < len(t) - lead and t[len(t) - 1 - trail] == "-":
    trail += 1
  extra_5, extra_3 = lead, trail
  subst = ins = dele = amb = n_bases = 0
  site_errors = 0
  tpos = start
  for i in range(lead, len(q) - trail):
    a, b = q[i], t[i]
    in_site = site is not None and site[0] <= tpos < site[1]
    if a == "-":
      dele += 1
      site_errors += in_site
      tpos += 1
    elif b == "-":
      ins += 1
      site_errors += in_site
    else:
      if a == "N":
        n_bases += 1
      elif a not in "ACGT":
        amb += 1
      elif a != b:
        subst += 1
        site_errors += in_site
      tpos += 1
  return {"cons_len": len(cons), "truth_len": len(truth),
          "missing_5": start, "missing_3": len(truth) - 1 - end,
          "extra_5": extra_5, "extra_3": extra_3, "subst": subst,
          "ins": ins, "del": dele, "ambiguous": amb, "n_bases": n_bases,
          "site_errors": site_errors if site is not None else "NA"}


COLUMNS = ["cons_len", "truth_len", "missing_5", "missing_3", "extra_5",
           "extra_3", "subst", "ins", "del", "ambiguous", "site_errors",
           "n_bases"]

# the genes of the segments, in the names of the consensus files of IRMA
GENES = {1: "PB2", 2: "PB1", 3: "PA", 4: "HA", 5: "NP", 6: "NA", 7: "MP",
         8: "NS"}


def main():
  args = [a for a in sys.argv[1:] if a not in ("--no-header", "--plain")]
  plain = "--plain" in sys.argv
  if len(args) != 3:
    sys.exit(__doc__)
  simdir, irmadir, label = args
  truth = {}
  for header, seq in read_fasta(os.path.join(simdir, "genome.fasta")):
    truth[segment_number(header)] = (header.split()[0], seq)
  site = None
  with open(os.path.join(simdir, "site.tsv")) as f:
    f.readline()
    fields = f.readline().split("\t")
    site_segment, site_start = fields[0], int(fields[1])
    site_len = len(fields[2])
  if "--no-header" not in sys.argv:
    print("\t".join(["label", "segment"] + COLUMNS))
  for n in sorted(truth):
    name, tseq = truth[n]
    if plain:
      files = glob.glob(os.path.join(irmadir, "A_%s*.fasta" % GENES[n]))
    else:
      files = glob.glob(os.path.join(irmadir, "amended_consensus",
                                     "*_%d.fa" % n))
    if not files:
      print("\t".join([label, str(n)] + ["NA"] * len(COLUMNS)))
      continue
    cons = read_fasta(files[0])[0][1].replace("-", "")
    s = (site_start, site_start + site_len) if name == site_segment \
      else None
    m = compare(tseq, cons, s)
    print("\t".join([label, str(n)] + [str(m[c]) for c in COLUMNS]))


if __name__ == "__main__":
  main()
