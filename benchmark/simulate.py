# (c) 2026 Giorgio Gonnella
"""Simulate reads of the influenza A multisegment RT-PCR, with a truth.

The reads are simulated from a complete genome (the eight segments with
their conserved termini; see benchmark/README.md for the one used), as
amplicons of the default scheme of flucut: the degenerate untailed Uni12
primer (AGCRAAAGCAGG) at the start of each segment, the tailed Uni13
primer (mbtuni-13) at its end. The bases under a primer come from the
primer, not from the segment: in particular the R of Uni12 is A or G at
random in each read, whatever the segment has.

Internal priming: in the HA segment (the record with "segment 4" in its
header), the position most similar to the Uni12 primer (fewest
mismatches, without gaps; found in the sequence, not chosen) is used as
an internal priming site: a share of the HA reads start there
with the bases of the primer, which replace those of the template.

With --barcodes, each amplicon gets the native barcoding adapter of
Oxford Nanopore (flank, barcode, flank: sequences of Dorado) of one of
the first N barcodes at both ends, as in a ligation library. Each read
gets, independently: an orientation; random flanking bases
(standing for adapter and barcode remnants: random, not a real adapter
sequence); sometimes a start or an end inside a primer, or a break in
the middle; and sequencing errors (substitutions, insertions and
deletions) at a rate drawn for each read. These are modelling choices,
documented here and in benchmark/README.md, not measured properties of
a run.

Outputs, in <outdir>:
  reads.fastq.gz   the reads
  truth.tsv.gz     for each read: name, segment, kind (full, internal,
                   dimer, broken), orientation, and the origin of each of
                   its bases (a string as long as the read: f flank,
                   including the barcoding adapters,
                   t primer tail, p primer bases at a terminus, i primer
                   bases at the internal site, T template), and the
                   barcode (or none)
  genome.fasta     the truth (a copy of the input genome)
  barcodes.fasta   with --barcodes: the barcoding adapters used
                   (flank, barcode, flank), named NB01 to NB04
  site.tsv         the internal priming site (segment, 0-based start,
                   template bases replaced, primer, mismatches)

Usage:
  python3 -I benchmark/simulate.py <genome.fasta> <outdir> \\
      [--reads-per-segment N] [--internal-share X] [--seed S]
"""

import argparse
import gzip
import os
import random

UNI12 = "AGCRAAAGCAGG"
MBT_TAIL = "ACGCGTGATC"
UNI13 = "AGTAGAAACAAGG"
# the native barcoding adapters of Oxford Nanopore (flank, barcode,
# flank), as defined in Dorado (utils/barcode_kits.cpp, commit e72f149,
# quoted in src/kits/kits.rs of Barbell 0.3.3): NB01 to NB04
NB_FRONT = "ATTGCTAAGGTTAA"
NB_REAR = "CAGCACCT"
NB_BARCODES = ["CACAAAGACACCGACAACTTTCTT", "ACAGACGACTACAAACGGAATCGA",
               "CCTGGTAACTGGGACACAAGACTC", "TAGGGAAACACGATAGAATCCGAA"]

COMPLEMENT = str.maketrans("ACGTRN", "TGCAYN")
IUPAC = {"R": "AG", "N": "ACGT"}


def revcomp(seq):
  """Reverse complement (A, C, G, T, R, N)."""
  return seq.translate(COMPLEMENT)[::-1]


def read_fasta(path):
  """The records of a FASTA file, as (name, sequence, header) tuples."""
  records = []
  for line in open(path):
    line = line.strip()
    if line.startswith(">"):
      records.append([line[1:].split()[0], [], line[1:]])
    elif records:
      records[-1][1].append(line.upper())
  return [(n, "".join(s), h) for n, s, h in records]


def matches(p, b):
  """A primer base (possibly degenerate) matches a template base."""
  return b in IUPAC.get(p, p)


def best_site(primer, seq, skip):
  """The position of seq most similar to the primer, without gaps (the
  primer anneals with mismatches), outside the first and last <skip>
  bases (the termini): returns (mismatches, start). The best position
  must be unique."""
  n = len(primer)
  scores = []
  for i in range(skip, len(seq) - n - skip + 1):
    d = sum(not matches(p, b) for p, b in zip(primer, seq[i:i + n]))
    scores.append((d, i))
  best = min(scores)
  ties = [i for d, i in scores if d == best[0]]
  if len(ties) > 1:
    raise SystemExit("several equally good sites: %s" % ties)
  return best


def resolve(primer, rng):
  """A copy of a degenerate primer, as one molecule of the mix."""
  return "".join(rng.choice(IUPAC[c]) if c in IUPAC else c for c in primer)


def add_errors(seq, labels, rate, rng):
  """Substitutions, insertions and deletions at the given total rate
  (40%, 30%, 30%); an inserted base takes the label of the base after
  it."""
  out, out_labels = [], []
  for c, lab in zip(seq, labels):
    r = rng.random()
    if r < rate * 0.4:
      out.append(rng.choice([b for b in "ACGT" if b != c]))
      out_labels.append(lab)
    elif r < rate * 0.7:
      out.append(rng.choice("ACGT"))
      out_labels.append(lab)
      out.append(c)
      out_labels.append(lab)
    elif r < rate:
      continue
    else:
      out.append(c)
      out_labels.append(lab)
  return "".join(out), "".join(out_labels)


def main():
  parser = argparse.ArgumentParser(
      description="Simulate Uni12/Uni13 amplicon reads with a truth.")
  parser.add_argument("genome")
  parser.add_argument("outdir")
  parser.add_argument("--reads-per-segment", type=int, default=400)
  parser.add_argument(
      "--internal-share", type=float, default=0.6,
      help="share of the HA reads started at the internal site "
           "(default 0.6: enough for the artefact to reach a majority "
           "consensus when nothing is trimmed)")
  parser.add_argument("--seed", type=int, default=1)
  parser.add_argument(
      "--barcodes", type=int, default=0, choices=range(0, 5),
      help="native barcoding: each amplicon gets the adapter of one of the "
           "first N barcodes (NB01-NB04) at both ends (default 0: none)")
  args = parser.parse_args()
  rng = random.Random(args.seed)
  genome = read_fasta(args.genome)
  os.makedirs(args.outdir, exist_ok=True)

  # the segments, named by their number (segment 4 is HA)
  segments = [(name, seq) for name, seq, _ in genome]
  ha = [(n, s) for n, s, h in genome if "segment 4" in h]
  if len(ha) != 1:
    raise SystemExit("one HA record (\"segment 4\" in its header) is "
                     "needed in the genome")
  ha_name, ha_seq = ha[0]
  dist, site = best_site(UNI12, ha_seq, 50)
  with open(os.path.join(args.outdir, "site.tsv"), "w") as f:
    f.write("segment\tstart\treplaced\tprimer\tmismatches\n")
    f.write("%s\t%d\t%s\t%s\t%d\n"
            % (ha_name, site, ha_seq[site:site + len(UNI12)], UNI12, dist))
  with open(os.path.join(args.outdir, "genome.fasta"), "w") as f:
    for _, seq, header in genome:
      f.write(">%s\n%s\n" % (header, seq))

  if args.barcodes > 0:
    with open(os.path.join(args.outdir, "barcodes.fasta"), "w") as f:
      for b in range(args.barcodes):
        f.write(">NB%02d\n%s\n" % (b + 1,
                                    NB_FRONT + NB_BARCODES[b] + NB_REAR))
  rc13 = revcomp(MBT_TAIL + UNI13)
  rc13_labels = "p" * len(UNI13) + "t" * len(MBT_TAIL)
  fq = gzip.open(os.path.join(args.outdir, "reads.fastq.gz"), "wt")
  truth = gzip.open(os.path.join(args.outdir, "truth.tsv.gz"), "wt")
  truth.write("read\tsegment\tkind\tstrand\tlabels\tbarcode\n")
  n = 0
  for name, seq in segments:
    for k in range(args.reads_per_segment):
      p12 = resolve(UNI12, rng)
      kind = "full"
      r = rng.random()
      if name == ha_name and r < args.internal_share:
        kind = "internal"
        body = seq[site + len(UNI12):len(seq) - len(UNI13)]
        left_labels = "i" * len(p12)
      else:
        body = seq[len(UNI12):len(seq) - len(UNI13)]
        left_labels = "p" * len(p12)
        if rng.random() < 0.01:
          kind = "dimer"
          body = ""
      amp = p12 + body + rc13
      labels = left_labels + "T" * len(body) + rc13_labels
      # a start or an end inside a primer, or a break in the middle
      r = rng.random()
      if r < 0.1:
        cut = rng.randint(1, 8)
        amp, labels = amp[cut:], labels[cut:]
      elif r < 0.2:
        cut = rng.randint(1, 15)
        amp, labels = amp[:-cut], labels[:-cut]
      elif r < 0.25 and kind != "dimer":
        cut = rng.randint(len(amp) // 4, 3 * len(amp) // 4)
        if rng.random() < 0.5:
          amp, labels = amp[:cut], labels[:cut]
        else:
          amp, labels = amp[cut:], labels[cut:]
        kind = "broken" if kind == "full" else kind + "_broken"
      barcode = "none"
      if args.barcodes > 0:
        b = rng.randrange(args.barcodes)
        barcode = "NB%02d" % (b + 1)
        adapter = NB_FRONT + NB_BARCODES[b] + NB_REAR
        amp = adapter + amp + revcomp(adapter)
        labels = "f" * len(adapter) + labels + "f" * len(adapter)
      strand = "+"
      if rng.random() < 0.5:
        amp, labels = revcomp(amp), labels[::-1]
        strand = "-"
      lf, rf = rng.randint(0, 40), rng.randint(0, 40)
      amp = ("".join(rng.choice("ACGT") for _ in range(lf)) + amp +
             "".join(rng.choice("ACGT") for _ in range(rf)))
      labels = "f" * lf + labels + "f" * rf
      rate = min(max(rng.gauss(0.04, 0.02), 0.005), 0.15)
      read, read_labels = add_errors(amp, labels, rate, rng)
      qual = "".join(chr(33 + rng.randint(8, 30)) for _ in read)
      n += 1
      rname = "r%05d" % n
      fq.write("@%s\n%s\n+\n%s\n" % (rname, read, qual))
      truth.write("%s\t%s\t%s\t%s\t%s\t%s\n"
                  % (rname, name, kind, strand, read_labels, barcode))
  fq.close()
  truth.close()


if __name__ == "__main__":
  main()
