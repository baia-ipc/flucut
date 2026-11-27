# (c) 2026 Giorgio Gonnella
"""Tables of the mixed samples (Markdown), from the tables of
benchmark/mixed_samples.py.

From second_subtype.tsv, what flucut reports for the HA and the NA
(the reads of a second subtype, above its threshold or not):

- one_subtype: the samples with one subtype in the gene (one strain,
  or two strains of the same subtype): how many are above the
  threshold, and the largest numbers seen;
- remarked: for each pair of strains of different subtypes, the runs
  (error rates and seeds) above the threshold with the right second
  subtype, by minor fraction;
- by_fraction: the same, all pairs together.

From agreement.tsv, the statistic of each of the eight segments, by
reference: the excess of reads agreeing on a reference other than the
leading one (the reads with several voting probes which all vote for
it, minus those expected if the probes voted independently of their
read; the largest excess among the other references). A mix is
detected in a segment if its statistic is above the highest one of the
samples of one strain.

- single: the statistic in the samples of one strain, by segment;
- detected: for each pair of strains and segment, the runs in which
  the mix is detected, by minor fraction;
- lowest: for each pair and segment, the lowest minor fraction detected
  in all runs.

Usage:
  python3 benchmark/mixed_samples_tables.py <workdir> \\
      one_subtype|remarked|by_fraction|single|detected|lowest
"""

import collections
import csv
import sys

GENES = [None, "PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]


def load(workdir):
  """The statistic (excess reads, reads agreeing, reads with several
  voting probes) by (major, minor, fraction) and segment, a value for
  each run."""
  samples = collections.defaultdict(lambda: collections.defaultdict(list))
  for r in csv.DictReader(open(workdir + "/agreement.tsv"), delimiter="\t"):
    samples[(r["major"], r["minor"], float(r["minor_fraction"]),
             r["error_rate"], r["seed"])][int(r["segment"])].append(r)
  stats = collections.defaultdict(lambda: collections.defaultdict(list))
  for (major, minor, fraction, _, _), segments in samples.items():
    for s, refs in segments.items():
      lead = max(refs, key=lambda r: int(r["votes"]))
      excess, agreeing = max(
          (int(r["agreeing"]) - float(r["expected"]), int(r["agreeing"]))
          for r in refs if r is not lead)
      stats[(major, minor, fraction)][s].append(
          (excess, agreeing, int(refs[0]["reads"])))
  return stats


def table(header, rows):
  print("| " + " | ".join(header) + " |")
  print("|" + "---|" * len(header))
  for row in rows:
    print("| " + " | ".join(str(v) for v in row) + " |")
  print()


def thresholds(stats):
  """The highest statistic of the samples of one strain, by segment."""
  highest = collections.defaultdict(float)
  for (_, minor, _), segments in stats.items():
    if minor == "none":
      for s, values in segments.items():
        highest[s] = max([highest[s]] + [v[0] for v in values])
  return highest


def single_table(stats):
  rows = []
  for s in range(1, 9):
    values = [v for (_, minor, _), segments in stats.items()
              if minor == "none" for v in segments[s]]
    reads = sorted(v[2] for v in values)
    rows.append([GENES[s], len(values), reads[len(reads) // 2],
                 max(v[1] for v in values),
                 "%.1f" % max(v[0] for v in values)])
  table(["segment", "samples", "reads with several voting probes (median)",
         "most reads agreeing on another reference", "highest excess"],
        rows)


def mixes(stats):
  pairs = []
  for major, minor, _ in stats:
    if minor != "none" and (major, minor) not in pairs:
      pairs.append((major, minor))
  fractions = sorted({f for _, minor, f in stats if minor != "none"},
                     reverse=True)
  return pairs, fractions


def detected_table(stats):
  highest = thresholds(stats)
  pairs, fractions = mixes(stats)
  rows = []
  for major, minor in pairs:
    for f in fractions:
      row = [major, minor, "%g%%" % (100 * f)]
      for s in range(1, 9):
        values = stats[(major, minor, f)][s]
        row.append("%d/%d" % (sum(v[0] > highest[s] for v in values),
                              len(values)))
      rows.append(row)
  table(["major", "minor", "minor fraction"] + GENES[1:], rows)


def lowest_table(stats):
  highest = thresholds(stats)
  pairs, fractions = mixes(stats)
  rows = []
  for major, minor in pairs:
    row = [major, minor]
    for s in range(1, 9):
      lowest = "-"
      for f in fractions:
        values = stats[(major, minor, f)][s]
        if not values or not all(v[0] > highest[s] for v in values):
          break
        lowest = "%g%%" % (100 * f)
      row.append(lowest)
    rows.append(row)
  table(["major", "minor"] + GENES[1:], rows)


def load_subtypes(workdir):
  """The rows of second_subtype.tsv with counts."""
  with open(workdir + "/second_subtype.tsv") as f:
    return [r for r in csv.DictReader(f, delimiter="\t")
            if r["reads"] != "NA"]


def right(r):
  """Whether a row is above the threshold with the subtype of the minor
  strain as the second one (or as the leading one, in an even mix)."""
  return r["above_threshold"] == "1" and r["minor_subtype"] in (
      r["second"], r["leading"])


def one_subtype_table(rows):
  out = []
  for gene in SUBTYPE_GENES:
    values = [r for r in rows
              if r["gene"] == gene and r["minor_subtype"] == "none"]
    reads = sorted(int(r["reads"]) for r in values)
    out.append([gene, len(values), reads[len(reads) // 2],
                sum(r["above_threshold"] == "1" for r in values),
                max(int(r["second_reads"]) for r in values),
                "%.1f" % max(float(r["second_excess"]) for r in values)])
  table(["gene", "samples", "reads with several voting probes (median)",
         "above the threshold", "most reads of a second subtype",
         "highest excess"], out)


def remarked_table(rows, by_pair):
  """The runs above the threshold with the right second subtype, of
  the samples with two subtypes in the gene."""
  runs = collections.defaultdict(list)
  for r in rows:
    if r["minor_subtype"] != "none":
      pair = (r["major"], r["minor"]) if by_pair else ()
      runs[pair + ("%g%%" % (100 * float(r["minor_fraction"])),
                   r["gene"])].append(right(r))
  keys = []
  for key in runs:
    if key[:-1] not in keys:
      keys.append(key[:-1])
  out = []
  for key in keys:
    cells = [runs.get(key + (gene,)) for gene in SUBTYPE_GENES]
    out.append(list(key) + ["%d/%d" % (sum(c), len(c)) if c else "-"
                            for c in cells])
  table((["major", "minor"] if by_pair else []) +
        ["minor fraction"] + SUBTYPE_GENES, out)


SUBTYPE_GENES = ["HA", "NA"]
SUBTYPE_TABLES = {
    "one_subtype": one_subtype_table,
    "remarked": lambda rows: remarked_table(rows, True),
    "by_fraction": lambda rows: remarked_table(rows, False)}


def main():
  if len(sys.argv) != 3:
    sys.exit(__doc__)
  if sys.argv[2] in SUBTYPE_TABLES:
    SUBTYPE_TABLES[sys.argv[2]](load_subtypes(sys.argv[1]))
    return
  stats = load(sys.argv[1])
  {"single": single_table, "detected": detected_table,
   "lowest": lowest_table}[sys.argv[2]](stats)


if __name__ == "__main__":
  main()
