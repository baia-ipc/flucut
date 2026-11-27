# (c) 2026 Giorgio Gonnella
"""The origin of the read ends of flucut runs, from their JSON summaries.

For each summary (flucut_summary.json of a run with references), one
line with the reads, the 5' ends by origin (terminal, internal, not
located), the share of the internal ones among the located, the
reference chosen for the HA with the share of internally primed reads
flucut reports for that segment, and the largest internal priming site
(reference, end, positions, reads). With --sites <n>, the <n> largest
sites of each summary instead, one per line.

Usage:
  python3 -I benchmark/read_origins.py [--sites <n>] <label>=<summary>...
Output: a TSV table on the standard output.
"""

import json
import sys


def main():
  args = sys.argv[1:]
  n_sites = 0
  if args[:1] == ["--sites"]:
    n_sites = int(args[1])
    args = args[2:]
  if not args:
    sys.exit(__doc__)
  if n_sites:
    print("label\treference\tend\tfrom\tto\treads")
  else:
    print("label\treads_in\treads_written\tleft_terminal\tleft_internal\t"
          "left_not_located\tinternal_share\tsubtype\tha_reference\t"
          "ha_internal\tlargest_site\tlargest_site_reads")
  for arg in args:
    label, _, path = arg.partition("=")
    with open(path) as f:
      everything = json.load(f)["all"]
    sites = sorted(everything.get("internal_sites", []),
                   key=lambda s: -s["reads"])
    if n_sites:
      for s in sites[:n_sites]:
        print("%s\t%s\t%s\t%d\t%d\t%d" % (label, s["reference"], s["end"],
                                          s["from"], s["to"], s["reads"]))
      continue
    left = everything["ends"]["left"]
    located = left["terminal"] + left["internal"]
    share = "%.3f" % (left["internal"] / located) if located else "NA"
    sample = everything.get("sample", {})
    ha = [s for s in sample.get("segments", []) if s["gene"] == "HA"]
    site = "%s:%s:%d-%d" % (sites[0]["reference"], sites[0]["end"],
                            sites[0]["from"], sites[0]["to"]) \
        if sites else "NA"
    print("\t".join(str(x) for x in [
        label, everything["reads"]["in"], everything["reads"]["written"],
        left["terminal"], left["internal"], left["unknown"], share,
        sample.get("subtype", "NA"),
        ha[0]["reference"] if ha else "NA",
        ha[0]["internal"] if ha else "NA", site,
        sites[0]["reads"] if sites else 0]))


if __name__ == "__main__":
  main()
