# (c) 2026 Giorgio Gonnella
"""Build a reference set for flucut from the seed sequences of IRMA.

IRMA assembles the reads against seed sequences which hold the coding
sequence of each segment only. flucut needs references which begin and
end with the conserved termini of the segments, where the Uni12 and
Uni13 primers bind: a read end at a segment terminus is recognised by
mapping within --terminal-window bases of a reference end, so a
reference missing its untranslated regions (UTRs) would turn genuine
segment termini into internal priming.

This script joins the two: to each seed it attaches the UTRs of a
complete record of the same segment and, for the HA and the NA, of the
same subtype. The result is a reference whose coding sequence is
IRMA's, so that both tools work on the same sequences, and whose ends
reach the termini, so that flucut can place the read ends.

The donors are by default the built-in reference set of flucut, which
holds a complete record of every HA subtype (H1 to H16), of every NA
subtype (N1 to N9) and of the other six segments. The UTRs of a
segment are nearly the same in every influenza A virus, but their
length differs between the subtypes of the HA and of the NA (the 5'
UTR of the HA is 17 bases in H14 and 43 in H2), which is why the donor
of an HA or NA seed must be of its subtype: a donor of another subtype
would put the start of the coding sequence at the wrong distance from
the terminus.

A reference built this way is not a record of a database: its coding
sequence and its UTRs come from two different viruses. Its description
says so, and the provenance table lists both sources. Use it to trim
the reads, never as a reference for an assembly.

The paths can be stored in the configuration file of flucut
(--save-config): a run without -r then uses this set, and the set can
be built again with --from-config when IRMA updates its seeds.

Usage:
  python3 scripts/irma_references.py [options]

Examples:
  # build the set from an IRMA installation and use it from now on
  python3 scripts/irma_references.py --irma /opt/flu-amd/IRMA_RES \\
      --save-config

  # build it again, with the paths stored before (after an IRMA update)
  python3 scripts/irma_references.py --from-config
"""

import argparse
import os
import re
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])

from fetch_references import log  # noqa: E402

# the conserved termini of the segments: a reference must begin with the
# first and end with the second (the fourth base of the 5' terminus is A
# or G, hence one mismatch allowed)
TERMINUS_5 = "AGCAAAAGCAGG"
TERMINUS_3 = "CCTTGTTTCTACT"
ONE_OFF = 1

# the segment of each gene name of a seed (A_HA_H5, A_MP, A_NS)
SEGMENT_OF = {"PB2": 1, "PB1": 2, "PA": 3, "HA": 4, "NP": 5, "NA": 6,
              "MP": 7, "M": 7, "NS": 8}
GENE_OF = {1: "polymerase PB2 (PB2)", 2: "polymerase PB1 (PB1)",
           3: "polymerase PA (PA)", 4: "hemagglutinin (HA)",
           5: "nucleoprotein (NP)", 6: "neuraminidase (NA)",
           7: "matrix protein (MP)", 8: "nonstructural protein (NS)"}

# bases of the seed matched to the donor to find the coding sequence in
# it, the share of them allowed to differ (the donor is another virus:
# the coding sequence differs by substitutions and, rarely, by indels
# at its ends) and where along the seed the match is tried: the first
# bases of a coding sequence can differ too much (the signal peptide of
# the HA), and a base inserted or deleted among them shifts the match,
# so the next windows are tried in turn
ANCHOR = 30
ANCHOR_MAX_ERROR_RATE = 0.3
ANCHOR_OFFSETS = [0, 15, 30, 60, 120, 240]


def parse_args():
  p = argparse.ArgumentParser(
    description="Build references for flucut from the seeds of IRMA.")
  p.add_argument("--irma", default="",
                 help="IRMA installation or its IRMA_RES directory; the "
                      "seeds are read from "
                      "modules/<module>/reference/consensus.fasta")
  p.add_argument("--seeds", default="",
                 help="the seed sequences (FASTA), instead of --irma")
  p.add_argument("--module", default="FLU",
                 help="the IRMA module whose seeds are used "
                      "(default: FLU)")
  p.add_argument("--donors", default="",
                 help="complete references giving the UTRs (FASTA); "
                      "default: the built-in set of flucut, from its "
                      "cache")
  p.add_argument("--output", default="",
                 help="the reference set written (FASTA); default: "
                      "irma-<module>-references.fasta in the cache of "
                      "flucut")
  p.add_argument("--provenance", default="",
                 help="table of the sources joined (TSV); default: next "
                      "to the output")
  p.add_argument("--save-config", action="store_true",
                 help="store the paths in the configuration file of "
                      "flucut, so that runs without -r use this set")
  p.add_argument("--from-config", action="store_true",
                 help="take the paths stored by --save-config (to build "
                      "the set again after an IRMA update)")
  return p.parse_args()


def cache_dir():
  """The cache directory of flucut, as the program finds it."""
  xdg = os.environ.get("XDG_CACHE_HOME")
  if xdg:
    return os.path.join(xdg, "flucut")
  return os.path.join(os.path.expanduser("~"), ".cache", "flucut")


def config_path():
  """The configuration file of flucut, as the program finds it."""
  xdg = os.environ.get("XDG_CONFIG_HOME")
  base = xdg if xdg else os.path.join(os.path.expanduser("~"), ".config")
  return os.path.join(base, "flucut", "config")


def read_config():
  """The keys of the configuration file of flucut."""
  out = {}
  try:
    for line in open(config_path()):
      line = line.split("#")[0]
      if "=" in line:
        key, value = line.split("=", 1)
        if value.strip():
          out[key.strip()] = value.strip()
  except FileNotFoundError:
    pass
  return out


def write_config(values):
  """Store <values> in the configuration file, keeping the other keys."""
  config = read_config()
  config.update(values)
  os.makedirs(os.path.dirname(config_path()), exist_ok=True)
  with open(config_path(), "w") as handle:
    handle.write("# the configuration of flucut (see src/user-config.h)\n")
    for key, value in sorted(config.items()):
      handle.write(f"{key} = {value}\n")
  log(f"the paths are stored in {config_path()}")


def read_fasta(path):
  """The records of a FASTA file, as (name, description, sequence)."""
  out, name, description, parts = [], None, "", []
  for line in open(path):
    if line.startswith(">"):
      if name is not None:
        out.append((name, description, "".join(parts).upper()))
      fields = line[1:].rstrip("\n").split(None, 1)
      name = fields[0]
      description = fields[1] if len(fields) > 1 else ""
      parts = []
    elif name is not None:
      parts.append(line.strip())
  if name is not None:
    out.append((name, description, "".join(parts).upper()))
  return out


def seeds_path(args):
  """The seed sequences of the IRMA module."""
  if args.seeds:
    return args.seeds
  if not args.irma:
    sys.exit("give --irma or --seeds (or --from-config)")
  for base in (args.irma, os.path.join(args.irma, "IRMA_RES")):
    path = os.path.join(base, "modules", args.module, "reference",
                        "consensus.fasta")
    if os.path.exists(path):
      return path
  sys.exit(f"no seeds of the module {args.module} under {args.irma} "
           f"(expected modules/{args.module}/reference/consensus.fasta)")


def seed_segment(name):
  """The segment and the subtype of a seed name (A_HA_H5 is segment 4
  of subtype H5, A_MP segment 7 of no subtype); (0, "") for a name
  which is not a segment of influenza A (e.g. the seeds of influenza
  B)."""
  fields = name.split("_")
  if not fields or fields[0].upper() != "A":
    return (0, "")
  gene = fields[1].upper() if len(fields) > 1 else ""
  subtype = fields[2].upper() if len(fields) > 2 else ""
  return (SEGMENT_OF.get(gene, 0), subtype)


def donor_subtype(description):
  """The subtype in the description of a donor record."""
  found = re.search(r"H\d+N\d+", description)
  return found.group(0) if found else ""


def donor_segment(description):
  """The segment in the description of a donor record."""
  found = re.search(r"segment (\d)", description)
  return int(found.group(1)) if found else 0


def mismatches(text, motif):
  return sum(a != b for a, b in zip(text, motif)) + \
    abs(len(text) - len(motif))


def best_match(sequence, anchor, start, stop):
  """The position in <sequence[start:stop]> where <anchor> matches with
  fewest mismatches, and the mismatches there; (-1, len) if none is
  within ANCHOR_MAX_ERROR_RATE."""
  allowed = int(len(anchor) * ANCHOR_MAX_ERROR_RATE)
  best, fewest = -1, len(anchor)
  for i in range(max(start, 0), min(stop, len(sequence) - len(anchor)) + 1):
    errors = mismatches(sequence[i:i + len(anchor)], anchor)
    if errors < fewest:
      best, fewest = i, errors
      if errors == 0:
        break
  return (best, fewest) if fewest <= allowed else (-1, len(anchor))


def anchored(donor, seed, at_start, room):
  """Where the coding sequence of the seed begins (at_start) or ends in
  the donor, from the first window of ANCHOR_OFFSETS which matches;
  (-1, 0, 0) if none does."""
  for offset in ANCHOR_OFFSETS:
    if offset + ANCHOR > len(seed):
      break
    if at_start:
      anchor = seed[offset:offset + ANCHOR]
      found, errors = best_match(donor, anchor, 0, room + offset)
      if found >= 0:
        return (found - offset, errors, offset)
    else:
      anchor = seed[-(offset + ANCHOR):len(seed) - offset or None]
      found, errors = best_match(donor, anchor,
                                 len(donor) - room - offset - ANCHOR,
                                 len(donor))
      if found >= 0:
        return (found + ANCHOR + offset, errors, offset)
  return (-1, 0, 0)


def utrs_of(donor, seed):
  """The UTRs of the donor around the coding sequence of the seed, found
  by matching the ends of the seed in it; (None, None, reason) if they
  cannot be placed or the result does not reach both termini."""
  room = 300    # the UTRs are at most some tens of bases long
  head, head_errors, head_offset = anchored(donor, seed, True, room)
  if head < 0:
    return (None, None, "the start of the coding sequence is not in the "
            "donor")
  tail, tail_errors, tail_offset = anchored(donor, seed, False, room)
  if tail < 0:
    return (None, None, "the end of the coding sequence is not in the "
            "donor")
  if head < 0 or tail > len(donor) or head >= tail:
    return (None, None, "the coding sequence is placed the wrong way "
            "round in the donor")
  utr5, utr3 = donor[:head], donor[tail:]
  built = utr5 + seed + utr3
  if mismatches(built[:len(TERMINUS_5)], TERMINUS_5) > ONE_OFF:
    return (None, None, "the result does not begin with the conserved "
            "terminus")
  if mismatches(built[-len(TERMINUS_3):], TERMINUS_3) > ONE_OFF:
    return (None, None, "the result does not end with the conserved "
            "terminus")
  log(f"    donor: coding sequence at {head}-{tail} ({head_errors} and "
      f"{tail_errors} mismatches in {ANCHOR} bases, matched {head_offset} "
      f"and {tail_offset} bases inside), UTRs {len(utr5)} and "
      f"{len(utr3)}")
  return (utr5, utr3, "")


def donors_for(donors, segment, subtype):
  """The donor records of the segment whose subtype fits <subtype> (the
  H of an HA seed, the N of an NA seed), in the order to try them; all
  the records of the segment if the seed has no subtype."""
  same = [d for d in donors if donor_segment(d[1]) == segment]
  if not subtype:
    return same
  out = []
  for record in same:
    found = donor_subtype(record[1])
    if not found:
      continue
    h, n = found.split("N")
    if (subtype.startswith("H") and h == subtype) or \
       (subtype.startswith("N") and "N" + n == subtype):
      out.append(record)
  return out


def build(seeds, donors, output, provenance):
  """Write the references and the provenance table; the records built
  and the seeds left out."""
  rows, written, skipped = [], 0, []
  with open(output + ".part", "w") as out:
    for name, _, seed in seeds:
      segment, subtype = seed_segment(name)
      if segment == 0:
        skipped.append((name, "not a segment of influenza A"))
        continue
      log(f"  {name}: segment {segment}"
          + (f", subtype {subtype}" if subtype else ""))
      candidates = donors_for(donors, segment, subtype)
      if not candidates:
        skipped.append((name, f"no donor of segment {segment}"
                        + (f" and subtype {subtype}" if subtype else "")))
        continue
      donor, utr5, utr3, problems = None, None, None, []
      for record in candidates:
        utr5, utr3, problem = utrs_of(record[2], seed)
        if not problem:
          donor = record
          break
        problems.append(f"{record[0]}: {problem}")
      if donor is None:
        skipped.append((name, "; ".join(problems)))
        continue
      sequence = utr5 + seed + utr3
      full = donor_subtype(donor[1])
      out.write(f">{name} segment {segment} {GENE_OF[segment]}: coding "
                f"sequence of the seed {name} of IRMA, untranslated "
                f"regions of {donor[0]}"
                + (f" ({full})" if full else "") + "\n")
      for i in range(0, len(sequence), 70):
        out.write(sequence[i:i + 70] + "\n")
      rows.append((name, segment, subtype, donor[0], full, len(utr5),
                   len(utr3), len(seed), len(sequence)))
      written += 1
  os.replace(output + ".part", output)
  with open(provenance, "w") as handle:
    handle.write("reference\tsegment\tseed_subtype\tutr_donor\t"
                 "donor_subtype\tutr5\tutr3\tcds_length\tlength\n")
    for row in rows:
      handle.write("\t".join(str(v) for v in row) + "\n")
  return written, skipped


def main():
  args = parse_args()
  config = read_config()
  if args.from_config:
    args.irma = args.irma or config.get("irma", "")
    args.module = config.get("irma_module", args.module)
  seeds_file = seeds_path(args)
  donors_file = args.donors or os.environ.get(
    "FLUCUT_BUILTIN_REFERENCES", "") or built_in(cache_dir())
  output = args.output or os.path.join(
    cache_dir(), f"irma-{args.module.lower()}-references.fasta")
  provenance = args.provenance or re.sub(r"\.fasta$", ".tsv", output)
  os.makedirs(os.path.dirname(output), exist_ok=True)
  log(f"seeds: {seeds_file}")
  log(f"donors of the untranslated regions: {donors_file}")
  seeds = read_fasta(seeds_file)
  donors = read_fasta(donors_file)
  written, skipped = build(seeds, donors, output, provenance)
  for name, why in skipped:
    log(f"left out: {name} ({why})")
  log(f"{written} references of the {len(seeds)} seeds written to "
      f"{output}")
  log(f"provenance written to {provenance}")
  if args.save_config:
    values = {"references": output}
    if args.irma:
      values["irma"] = args.irma
      values["irma_module"] = args.module
    write_config(values)
  if not written:
    sys.exit("no reference could be built")


def built_in(cache):
  """The built-in reference set in the cache of flucut, the newest
  version there."""
  found = sorted(f for f in os.listdir(cache)
                 if re.fullmatch(r"flu-a-references-v\d+\.fasta", f)) \
    if os.path.isdir(cache) else []
  if not found:
    sys.exit(f"no built-in reference set in {cache}: run flucut once to "
             f"download it, or give --donors")
  return os.path.join(cache, found[-1])


if __name__ == "__main__":
  main()
