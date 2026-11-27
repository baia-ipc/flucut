# (c) 2026 Giorgio Gonnella
"""Tables of the calibration of the reference choice (Markdown).

From accuracy.tsv and votes.tsv of benchmark/reference_choice.py:

- segments: for each simulated segment and reference, the vote share,
  the identity and the accuracy (the share of correct read ends), at
  one error rate;
- pairs: in the segments without a clear lead, each other reference
  compared with the leading one, by the ratio of their votes: what
  choosing it instead would cost;
- rule: what the rule of flucut chooses with other values of its three
  constants (benchmark/reference_rule.py), and the accuracy lost
  against the best reference of each segment; a segment whose
  reference is not chosen keeps all of them (the subtype "all" of
  accuracy.tsv);
- follow: in the segments where following the common subtype changes
  the reference, that reference against the leading one, by the ratio
  of their votes (the close share decides where the two differ);
- all: the reference chosen against keeping all the references of the
  segment;
- votes: how often a leading reference with few votes is the one which
  leads with all the reads;
- termini: the bases which each reference lacks at its ends, and the
  accuracy with it (its read ends of unknown origin): a reference
  which lacks more than the primer at an end cannot be told from an
  internal position there.

Usage:
  python3 benchmark/reference_choice_tables.py <workdir> \\
      <references.fasta> segments|pairs|rule|follow|all|votes|termini
"""

import collections
import csv
import re
import sys

import flucut_runs
import pairwise
import reference_rule
import score_read_ends

GENES = [None, "PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]
SUBTYPES = ["H1N1", "H3N2", "H5N1", "H7N9", "H9N2"]


def load(workdir):
  """accuracy[(genome, error rate, segment)][reference] = (identity,
  accuracy); votes[(genome, error rate, sample)][segment][reference];
  subtype[reference]."""
  accuracy = collections.defaultdict(dict)
  for r in csv.DictReader(open(workdir + "/accuracy.tsv"), delimiter="\t"):
    name = r["reference"] if r["subtype"] != "all" else "all:" + \
        r["reference"]
    accuracy[(r["genome"], r["error_rate"], int(r["segment"]))][name] = (
        float(r["identity"]), int(r["correct"]) / max(int(r["ends"]), 1))
  votes = collections.defaultdict(lambda: collections.defaultdict(dict))
  subtype = {}
  for r in csv.DictReader(open(workdir + "/votes.tsv"), delimiter="\t"):
    votes[(r["genome"], r["error_rate"], r["sample"])][
        int(r["segment"])][r["reference"]] = int(r["votes"])
    subtype[r["reference"]] = r["subtype"]
  return accuracy, votes, subtype


def table(header, rows):
  print("| " + " | ".join(header) + " |")
  print("|" + "---|" * len(header))
  for row in rows:
    print("| " + " | ".join(str(v) for v in row) + " |")
  print()


def pct(x):
  return "%.0f" % (100 * x)


def segments_table(accuracy, votes, subtype, error_rate="0.08"):
  """Vote share / identity / accuracy of each reference (percent); the
  reference chosen by flucut is in bold."""
  rows = []
  for (genome, e, sample), segs in sorted(votes.items()):
    if e != error_rate or sample != "even":
      continue
    chosen = reference_rule.decide(segs, subtype)
    for s in sorted(segs):
      total = sum(segs[s].values())
      cells = []
      for ref, n in segs[s].items():
        identity, acc = accuracy[(genome, e, s)][ref]
        cell = "%s / %s / %s" % (pct(n / max(total, 1)), pct(identity),
                                 pct(acc))
        cells.append("**%s**" % cell if ref == chosen[s][0] else cell)
      rows.append([genome, GENES[s], total] + cells + [chosen[s][1]])
  table(["genome", "segment", "votes"] + SUBTYPES + ["chosen"], rows)


def pairs_table(accuracy, votes, subtype, lacking):
  """Each other reference against the leading one, in the segments
  without a clear lead of the samples with 400 reads per segment."""
  bins = [(0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 1.01)]

  def whole(ref):
    return lacking[ref][0] <= 12 and lacking[ref][1] <= 13

  for title, keep in (
      ("all segments", lambda s, a, b: True),
      ("HA and NA", lambda s, a, b: s in (4, 6)),
      ("the other segments", lambda s, a, b: s not in (4, 6)),
      ("the other segments, both references whole up to the primers",
       lambda s, a, b: s not in (4, 6) and whole(a) and whole(b))):
    stats = {b: [] for b in bins}
    for (genome, e, sample), segs in votes.items():
      if sample != "even":
        continue
      for s, v in segs.items():
        best = reference_rule.leader(v)
        clear = reference_rule.choice(v, subtype, 20, 0.6)[2]
        if clear or v[best] < 20:
          continue
        for ref, n in v.items():
          if ref == best or n == 0 or not keep(s, best, ref):
            continue
          ratio = n / v[best]
          for b in bins:
            if b[0] <= ratio < b[1]:
              stats[b].append(accuracy[(genome, e, s)][ref][1] -
                              accuracy[(genome, e, s)][best][1])
    print("Segments without a clear lead, %s:\n" % title)
    rows = []
    for b in bins:
      d = sorted(stats[b])
      if not d:
        continue
      rows.append(["%.2f-%.2f" % (b[0], min(b[1], 1)), len(d),
                   "%+.1f" % (100 * sum(d) / len(d)),
                   "%+.1f" % (100 * d[len(d) // 2]),
                   sum(x < -0.02 for x in d), sum(x > 0.02 for x in d),
                   "%+.0f" % (100 * d[0])])
    table(["votes / votes of the leader", "pairs", "mean difference",
           "median", "worse by > 2", "better by > 2", "worst"], rows)


def regrets(accuracy, votes, subtype, constants, keep, all_refs):
  """The accuracy lost against the best reference of the segment, in
  the segments chosen as <how>, by how; a segment whose reference is
  not chosen keeps all of them."""
  lost = collections.defaultdict(list)
  for (genome, e, sample), segs in votes.items():
    for s, (ref, how) in reference_rule.decide(
        segs, subtype, **constants).items():
      if not keep(s):
        continue
      acc = accuracy[(genome, e, s)]
      best = max(a for ref, (_, a) in acc.items()
                  if not ref.startswith("all:"))
      lost[how].append(best - (acc[all_refs[(genome, e, s)]][1]
                               if ref is None else acc[ref][1]))
  return lost


def rule_table(accuracy, votes, subtype, all_refs):
  """The accuracy lost with other values of the constants."""
  variants = [("close share", "close_share", [None, 0.75, 0.5, 0.25, 0.0]),
              ("least votes", "min_votes", [5, 10, 20, 40, 80]),
              ("clear share", "min_share", [0.5, 0.6, 0.7, 0.8, 0.9])]
  for name, key, values in variants:
    print("%s:\n" % name)
    rows = []
    for value in values:
      constants = dict(reference_rule.DEFAULTS, **{key: value})
      lost = regrets(accuracy, votes, subtype, constants, lambda s: True,
                     all_refs)
      every = [x for how in lost for x in lost[how]]
      unclear = [x for how in ("follows", "leader") for x in lost[how]]
      row = ["never" if value is None else value,
             len(lost["clear"]), len(lost["follows"]), len(lost["leader"]),
             len(lost["none"])]
      for group in (unclear, lost["follows"], every):
        row.append("%.1f" % (100 * sum(group) / max(len(group), 1)))
      row.append(sum(x > 0.05 for x in every))
      rows.append(row)
    table([name, "clear", "follow", "leader", "all kept",
           "lost, not clear", "lost, follow", "lost, every segment",
           "segments losing > 5"], rows)


def follow_table(accuracy, votes, subtype, all_refs):
  """Where following the common subtype changes the reference, the
  accuracy of the followed reference against that of the leader."""
  bins = [(0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 1.01)]
  stats = {b: [] for b in bins}
  for (genome, e, sample), segs in votes.items():
    follows = reference_rule.decide(segs, subtype, close_share=0.0)
    for s, (ref, how) in follows.items():
      best = reference_rule.leader(segs[s])
      if how != "follows" or ref == best or segs[s][best] < 20:
        continue
      acc = accuracy[(genome, e, s)]
      ratio = segs[s][ref] / segs[s][best]
      for b in bins:
        if b[0] <= ratio < b[1]:
          stats[b].append(acc[ref][1] - acc[best][1])
  rows = []
  for b in bins:
    d = sorted(stats[b])
    if not d:
      continue
    rows.append(["%.2f-%.2f" % (b[0], min(b[1], 1)), len(d),
                 "%+.1f" % (100 * sum(d) / len(d)),
                 "%+.1f" % (100 * d[len(d) // 2]),
                 sum(x > 0.02 for x in d), sum(x < -0.02 for x in d),
                 "%+.0f" % (100 * d[0])])
  table(["votes / votes of the leader", "segments", "mean difference",
         "median", "better by > 2", "worse by > 2", "worst"], rows)


def all_table(accuracy, votes, subtype, all_refs):
  """The reference chosen against keeping all the references of the
  segment, by how it was chosen."""
  rows = collections.defaultdict(list)
  for (genome, e, sample), segs in votes.items():
    for s, (ref, how) in reference_rule.decide(segs, subtype).items():
      if ref is None:
        continue
      acc = accuracy[(genome, e, s)]
      rows[how].append(acc[ref][1] - acc[all_refs[(genome, e, s)]][1])
      rows["every segment"].append(rows[how][-1])
  table(["chosen", "segments", "mean difference to all the references",
         "better by > 2", "worse by > 2", "worst"],
        [[how, len(d), "%+.1f" % (100 * sum(d) / len(d)),
          sum(x > 0.02 for x in d), sum(x < -0.02 for x in d),
          "%+.0f" % (100 * min(d))]
         for how, d in sorted(rows.items())])


def votes_table(accuracy, votes, subtype):
  """Whether a leading reference with few votes is the one leading
  with all the reads, by its votes and share."""
  bins = [(1, 5), (5, 10), (10, 20), (20, 40), (40, 80), (80, 10**9)]
  shares = [(0.0, 0.5), (0.5, 0.6), (0.6, 0.7), (0.7, 0.8), (0.8, 1.01)]
  count = collections.defaultdict(lambda: [0, 0, 0.0])
  for (genome, e, sample), segs in votes.items():
    if sample == "even":
      continue
    full = votes[(genome, e, "even")]
    for s, v in segs.items():
      best = reference_rule.leader(v)
      if v[best] == 0:
        continue
      share = v[best] / sum(v.values())
      acc = accuracy[(genome, e, s)]
      for b in bins:
        for sh in shares:
          if b[0] <= v[best] < b[1] and sh[0] <= share < sh[1]:
            c = count[(b, sh)]
            c[0] += 1
            c[1] += best == reference_rule.leader(full[s])
            c[2] += max(a for _, a in acc.values()) - acc[best][1]
  rows = []
  for b in bins:
    row = ["%d-%d" % (b[0], b[1] - 1) if b[1] < 10**9 else ">= %d" % b[0]]
    for sh in shares:
      n, same, lost = count[(b, sh)]
      row.append("%d / %d (%.1f)" % (same, n, 100 * lost / n) if n else "")
    rows.append(row)
  table(["votes of the leader"] +
        ["share %d-%d%%" % (100 * a, min(100 * b, 100)) for a, b in shares],
        rows)


def missing(whole, seq):
  """The bases which seq lacks at its start, against a sequence of the
  same segment with its terminus: the shift at which the first bases
  of seq are most similar to those of the whole one."""
  def differences(shift):
    return sum(a != b for a, b in zip(whole[shift:shift + 40], seq[:40]))
  return min(range(41), key=differences)


def reference_termini(references):
  """The bases each reference lacks at its start and at its end, by
  reference name: against the first reference of its segment with both
  conserved termini."""
  refs = score_read_ends.read_fasta(references)
  lacking = {}
  for s in range(1, 9):
    same = [(name, seq) for name, seq, n in refs if n == s]
    whole = [seq for _, seq in same if seq.startswith("AGC") and
             seq.endswith("CTACT")][0]
    for name, seq in same:
      lacking[name] = (missing(whole, seq), missing(whole[::-1], seq[::-1]))
  return lacking


def termini_table(references, accuracy=None):
  """The bases each reference lacks at its ends (start, end), with the
  accuracy with it and its ends of unknown origin, in percent."""
  lacking = reference_termini(references)
  subtype, segment = {}, {}
  for line in open(references):
    if line.startswith(">"):
      name = line[1:].split()[0]
      subtype[name] = flucut_runs.SUBTYPE.search(line).group(0)
      segment[name] = int(re.search(r"segment (\d)", line).group(1))
  measured = collections.defaultdict(list)
  for (_, _, s), refs in (accuracy or {}).items():
    for ref, (_, a) in refs.items():
      if not ref.startswith("all:"):
        measured[ref].append(a)
  rows = []
  for s in range(1, 9):
    cells = {}
    for name in lacking:
      if segment[name] != s:
        continue
      cells[subtype[name]] = "%d, %d" % lacking[name]
      if measured.get(name):
        cells[subtype[name]] += " (%.0f)" % (
            100 * sum(measured[name]) / len(measured[name]))
    rows.append([GENES[s]] + [cells[st] for st in SUBTYPES])
  table(["segment"] + SUBTYPES, rows)


def termini_accuracy(references, accuracy):
  """The accuracy by the bases the reference lacks at its ends, in the
  internal genes (in HA and NA a reference of another subtype is too
  far for the probes, whatever its termini)."""
  lacking = reference_termini(references)
  bins = [(0, 0), (1, 12), (13, 10 ** 9)]
  groups = collections.defaultdict(list)
  for (_, _, s), refs in accuracy.items():
    if s in (4, 6):
      continue
    for ref, (_, a) in refs.items():
      if ref.startswith("all:"):
        continue
      for b in bins:
        if b[0] <= max(lacking[ref]) <= b[1]:
          groups[b].append(a)
  print("The internal genes (not HA and NA):\n")
  table(["bases lacking at an end", "references x samples",
         "mean accuracy"],
        [["%d-%d" % b if b[1] < 10 ** 9 else "> %d" % (b[0] - 1),
          len(groups[b]), "%.0f" % (100 * sum(groups[b]) / len(groups[b]))]
         for b in bins if groups[b]])


def main():
  if len(sys.argv) != 4:
    sys.exit(__doc__)
  workdir, references, what = sys.argv[1:]
  accuracy, votes, subtype = load(workdir)
  if what == "termini":
    termini_table(references, accuracy)
    termini_accuracy(references, accuracy)
    return
  all_refs = {}
  for key, refs in accuracy.items():
    for ref in refs:
      if ref.startswith("all:"):
        all_refs[key] = ref
  if what in ("rule", "all", "follow"):
    {"rule": rule_table, "all": all_table,
     "follow": follow_table}[what](accuracy, votes, subtype, all_refs)
    return
  if what == "pairs":
    pairs_table(accuracy, votes, subtype, reference_termini(references))
    return
  {"segments": segments_table,
   "votes": votes_table}[what](accuracy, votes, subtype)


if __name__ == "__main__":
  main()
