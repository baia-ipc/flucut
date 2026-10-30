# (c) 2026 Giorgio Gonnella
"""Count the codon of the Uni12-like HA site in raw Nanopore reads.

Reads are matched (both strands, exact) on a flank of the template
downstream of the 12-base site; the 12 bases before the flank are the
site, and its bases 5-7 are the codon (ATA template, AAA artefact).
Reads are classified by where the site lies in the read: near the read
start (<= --window bases: possible internal priming) or further in
(the read spans the site: the template).

Usage: count_codon_nanopore.py <isolate.fasta> <HA id> <site 1-based> <fastq.gz>...
"""
import collections, gzip, sys

def rc(s):
  return s[::-1].translate(str.maketrans("ACGT", "TGCA"))

def main():
  fasta, ha_id, site = sys.argv[1], sys.argv[2], int(sys.argv[3])
  seq, name = {}, None
  for l in open(fasta):
    l = l.strip()
    if l.startswith(">"): name = l[1:].split()[0]; seq[name] = []
    elif name: seq[name].append(l.upper())
  ha = "".join(seq[ha_id])
  s0 = site - 1
  flank = ha[s0 + 12:s0 + 24]
  # note: the reference may itself carry the artefact codon
  print("# reference site", ha[s0:s0 + 12], "codon", ha[s0 + 4:s0 + 7],
        "flank", flank)
  print("run\tclass\tcodon\treads")
  for path in sys.argv[4:]:
    c = collections.Counter()
    # a truncated download (partial gzip) is read up to its end
    with gzip.open(path, "rt") as f:
     try:
      for i, l in enumerate(f):
        if i % 4 != 1: continue
        r = l.strip().upper()
        for strand, q in (("+", r), ("-", rc(r))):
          j = q.find(flank)
          if j >= 12:
            cod = q[j - 12 + 4:j - 12 + 7]
            # internal priming starts the read at the site
            start = j - 12
            c[(strand, "site_within_60_of_read_start" if start <= 60 else "site_inside", cod)] += 1
     except EOFError:
      pass
    for (st, cl, cod), n in sorted(c.items(), key=lambda x: (x[0][0], x[0][1], -x[1])):
      print("%s\t%s%s\t%s\t%d" % (path.split("/")[-1][:12], st, " " + cl, cod, n))

main()
