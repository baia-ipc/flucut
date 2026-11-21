# (c) 2026 Giorgio Gonnella
"""Download the genomes of a table of accessions from NCBI.

The table (e.g. benchmark/reference_choice_genomes.tsv) has a line for
each genome: label, strain, subtype, the accessions of its 8 segments
(comma-separated, segments 1 to 8), and why it was chosen. Each genome
is written to <outdir>/<label>.fasta, with the segments in order and
"segment <n>" in each header (as benchmark/simulate.py needs). A genome
is refused if a record does not start with the Uni12 terminus and end
with the Uni13 terminus (at most one mismatch each): the simulation
needs complete segments.

The requests carry no name or email, only a generic client name.

Usage:
  python3 -I benchmark/fetch_genomes.py <table.tsv> <outdir>
"""

import os
import sys
import time
import urllib.parse
import urllib.request

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
STARTS = ["AGCAAAAGCAGG", "AGCGAAAGCAGG"]
END = "CCTTGTTTCTACT"


def mismatches(a, b):
  return sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))


def read_table(path):
  """The genomes of the table, as (label, strain, subtype, accessions)."""
  genomes = []
  for line in open(path):
    if line.startswith("#") or line.startswith("label\t"):
      continue
    f = line.rstrip("\n").split("\t")
    if len(f) >= 4:
      genomes.append((f[0], f[1], f[2], f[3].split(",")))
  return genomes


def fetch(accessions):
  """The sequences of the accessions, in their order."""
  data = urllib.parse.urlencode({
      "db": "nuccore", "rettype": "fasta", "retmode": "text",
      "id": ",".join(accessions)}).encode()
  request = urllib.request.Request(
      EFETCH, data=data, headers={"User-Agent": "python-urllib"})
  for wait in (1, 4, 16, None):
    try:
      text = urllib.request.urlopen(request, timeout=120).read().decode()
      break
    except OSError:
      if wait is None:
        raise
      time.sleep(wait)
  records = {}
  for record in text.split(">")[1:]:
    header, _, body = record.partition("\n")
    records[header.split()[0]] = body.replace("\n", "").upper()
  missing = [a for a in accessions if a not in records]
  if missing:
    sys.exit("not returned by NCBI: %s" % ", ".join(missing))
  return [records[a] for a in accessions]


def complete(seq):
  """Whether a segment has both conserved termini."""
  return (min(mismatches(seq[:12], s) for s in STARTS) <= 1 and
          mismatches(seq[-13:], END) <= 1)


def main():
  if len(sys.argv) != 3:
    sys.exit(__doc__)
  table, outdir = sys.argv[1:]
  os.makedirs(outdir, exist_ok=True)
  for label, strain, subtype, accessions in read_table(table):
    path = os.path.join(outdir, label + ".fasta")
    if os.path.exists(path):
      continue
    if len(accessions) != 8:
      sys.exit("%s: 8 accessions are needed" % label)
    seqs = fetch(accessions)
    for n, (accession, seq) in enumerate(zip(accessions, seqs), 1):
      if not complete(seq):
        sys.exit("%s: segment %d (%s) lacks a conserved terminus"
                 % (label, n, accession))
    with open(path, "w") as out:
      for n, (accession, seq) in enumerate(zip(accessions, seqs), 1):
        out.write(">%s Influenza A virus (%s(%s)) segment %d\n%s\n"
                  % (accession, strain, subtype, n, seq))
    sys.stderr.write("%s: %s (%s)\n" % (label, strain, subtype))
    time.sleep(0.5)


if __name__ == "__main__":
  main()
