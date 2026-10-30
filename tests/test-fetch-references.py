# (c) 2026 Giorgio Gonnella
"""Tests of scripts/fetch_references.py.

Without arguments, checks offline the esearch terms built for the
eight influenza A segments: each segment is queried by gene names
which NCBI uses (segment 7 by M1/M2, segment 8 by NS1/NEP/NS2, as
MP[Gene] and NS[Gene] match no RefSeq record), and the --subtype and
--refseq-only terms are kept.

With --live, also downloads the RefSeq references of all eight
segments from NCBI (network access required, through HTTPS_PROXY if
set) and checks that every segment is in the FASTA output and in the
provenance table, with one row per sequence.

Usage (from the repository root):

  python3 -I tests/test-fetch-references.py           # offline
  python3 -I tests/test-fetch-references.py --live    # also NCBI
"""

import importlib.util
import os
import sys
import tempfile

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "..", "scripts", "fetch_references.py")

SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]


def load_script():
  """Import the script as a module (from its path: run with -I)."""
  spec = importlib.util.spec_from_file_location("fetch_references", SCRIPT)
  module = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(module)
  return module


def check(cond, msg):
  """Exit with status 1 and the message if the condition is false."""
  if not cond:
    sys.stderr.write("[FAIL] fetch_references: %s\n" % msg)
    sys.exit(1)


def test_queries(fr):
  """The esearch terms of the eight segments."""
  org = "Influenza A virus[Organism] AND "
  expected = {
    "PB2": "PB2[Gene]", "PB1": "PB1[Gene]", "PA": "PA[Gene]",
    "HA": "HA[Gene]", "NP": "NP[Gene]", "NA": "NA[Gene]",
    "MP": "(M1[Gene] OR M2[Gene])",
    "NS": "(NS1[Gene] OR NEP[Gene] OR NS2[Gene])",
  }
  check(fr.resolve_segments("all") == SEGMENTS,
        "--segment all is not the eight segments in order")
  for seg in fr.resolve_segments("all"):
    term = fr.build_query(seg, None, False)
    check(term == org + expected[seg], "%s: term %r" % (seg, term))
    check("MP[Gene]" not in term and "NS[Gene]" not in term,
          "%s: term %r uses a segment name as gene" % (seg, term))
  term = fr.build_query("MP", "H5N1", True)
  check(term == org + expected["MP"] + " AND H5N1 AND refseq[filter]",
        "MP with --subtype and --refseq-only: term %r" % term)
  check(fr.resolve_segments("ns") == ["NS"], "lowercase segment name")


def test_live(fr):
  """Download the RefSeq references of all segments from NCBI."""
  # the records downloaded for each segment
  found = {}
  fetch_segment = fr.fetch_segment
  def fetch_and_count(segment, *args):
    recs = fetch_segment(segment, *args)
    found[segment] = [r["accession"] for r in recs]
    return recs
  fr.fetch_segment = fetch_and_count
  with tempfile.TemporaryDirectory() as tmp:
    out = os.path.join(tmp, "all.fasta")
    status = fr.main(["--segment", "all", "--max", "4", "--refseq-only",
                      "--output", out])
    check(status == 0, "live run: exit status %d" % status)
    with open(out) as f:
      headers = [l[1:] for l in f if l.startswith(">")]
    with open(out + ".provenance.tsv") as f:
      rows = [l.rstrip("\n").split("\t") for l in f]
  check(rows[0] == ["accession", "description", "length", "source_url",
                    "fetch_date"], "provenance header %r" % rows[0])
  rows = rows[1:]
  check(len(rows) == len(headers),
        "%d provenance rows for %d sequences" % (len(rows), len(headers)))
  for (acc, desc, length, url, _), header in zip(rows, headers):
    check(header == "%s %s\n" % (acc, desc),
          "provenance row %s does not match the FASTA" % acc)
    check(int(length) > 0 and url.startswith(fr.EUTILS_BASE),
          "provenance row %s: length %s, url %s" % (acc, length, url))
  for num, seg in enumerate(SEGMENTS, 1):
    check(found.get(seg), "live run: no record of segment %d (%s)"
          % (num, seg))
    sys.stderr.write("segment %d (%s): %s\n"
                     % (num, seg, ", ".join(found[seg])))
  check(sum(len(a) for a in found.values()) == len(rows),
        "the provenance rows are not the records of the segments")


def main():
  fr = load_script()
  test_queries(fr)
  if "--live" in sys.argv[1:]:
    test_live(fr)
    print("[OK] fetch_references (with the NCBI download)")
  else:
    print("[OK] fetch_references")


if __name__ == "__main__":
  main()
