# (c) 2026 Giorgio Gonnella
"""Simulate reads of the influenza multisegment RT-PCR, with a truth.

The reads are simulated from a complete genome (the eight segments with
their conserved termini; see benchmark/README.md for the one used), as
amplicons of the default scheme of flucut for the type of the genome
(--virus).

For influenza A (the default) both primers are the same in every
segment: the degenerate untailed Uni12 primer (AGCRAAAGCAGG) at the
start of each segment and the tailed Uni13 primer (mbtuni-13) at its
end. For influenza B the primers are the IBV-GA2 cocktail, whose six
pairs each serve one or two segments, so the pair of a segment is the
one designed for it; each of its primers is a 5' tail and then the
bases of the segment at that end.

The bases under a primer come from the primer, not from the segment: in
particular the R of Uni12 is A or G at random in each read, whatever
the segment has.

Internal priming: in the HA segment (the record with "segment 4" in its
header; "segment 4" is the HA of both types), the position most similar
to the primer of the left end (fewest mismatches, without gaps; found
in the sequence, not chosen) is used as an internal priming site: a
share of the HA reads start there with the bases of the primer, which
replace those of the template.

With --barcodes, each amplicon gets the native barcoding adapter of
Oxford Nanopore (flank, barcode, flank: sequences of Dorado) of one of
the first N barcodes at both ends, as in a ligation library. Each read
gets, independently: an orientation; random flanking bases
(standing for adapter and barcode remnants: random, not a real adapter
sequence); sometimes a start or an end inside a primer, or a break in
the middle; and sequencing errors (substitutions, insertions and
deletions) at a rate drawn for each read (--error-rate: its mean).
With --joined-share, a share of the reads is two amplicons joined end to
end, as a ligation library also gives: the second one is of any segment,
each has its own orientation, and 0 to 20 other bases lie between them.
These are modelling choices, documented here and in
benchmark/README.md, not measured properties of a run.

Outputs, in <outdir>:
  reads.fastq.gz   the reads
  truth.tsv.gz     for each read: name, segment, kind (full, internal,
                   dimer, broken, joined), orientation, the origin of each of
                   its bases (a string as long as the read: f flank,
                   including the barcoding adapters,
                   t primer tail, p primer bases at a terminus which
                   are the sequence of the segment there, i primer
                   bases which are not (at the internal priming site,
                   or of a primer of other segments, see
                   --wrong-pair-share), T template), the
                   barcode (or none), and the number of amplicons in
                   the read (2 if joined)
  genome.fasta     the truth (a copy of the input genome)
  barcodes.fasta   with --barcodes: the barcoding adapters used
                   (flank, barcode, flank), named NB01 to NB04
  site.tsv         the internal priming site (segment, 0-based start,
                   template bases replaced, primer, mismatches)

The reads are written segment after segment; --shuffle writes them in
random order, as a sequencer does (flucut chooses its references from
the first reads). --segment-reads gives some segments another number of
reads (e.g. 1:8,2:20: 8 reads of segment 1, 20 of segment 2), for
samples with uneven coverage. --first-site accepts an HA with several
positions equally similar to the primer. --read-prefix names the reads
<prefix>r00001, ..., so that the reads of several simulations can be
put together (benchmark/mix_reads.py).

Usage:
  python3 -I benchmark/simulate.py <genome.fasta> <outdir> \\
      [--reads-per-segment N] [--internal-share X] [--seed S]
"""

import argparse
import gzip
import os
import random
import re

UNI12 = "AGCRAAAGCAGG"
MBT_TAIL = "ACGCGTGATC"
UNI13 = "AGTAGAAACAAGG"

# The primers of each segment, by type of influenza: the primer of the
# left end of an amplicon and the primer of its right end, each as it is
# ordered (5' -> 3'), with the length of its 5' tail (the rest is the
# sequence of the segment at that end). They are the primers built into
# flucut, see src/primers.h.
#
# Influenza A has one pair for every segment; the pairs of influenza B
# (the IBV-GA2 cocktail) each serve one or two segments, and the
# segments are numbered in the numbering of their type (influenza B
# numbers PB1 1 and PB2 2). The M segment of influenza B has two
# primers of the left end in the cocktail; the first is used here.
PRIMERS = {
  "a": {n: ((UNI12, 0), (MBT_TAIL + UNI13, len(MBT_TAIL)))
        for n in range(1, 9)},
  "b": {
      1: (("GGGGGGAGCAGAAGCGGAGC", 6),
          ("CCGGGTTATTAGTAGAAACACGAGC", 10)),       # PB1
      2: (("GGGGGGAGCAGAAGCGGAGC", 6),
          ("CCGGGTTATTAGTAGAAACACGAGC", 10)),       # PB2
      3: (("GGGGGGAGCAGAAGCGGTGC", 6),
          ("CCGGGTTATTAGTAGAAACACGTGC", 10)),       # PA
      4: (("GGGGGGAGCAGAAGCAGAGC", 6),
          ("CCGGGTTATTAGTAGTAACAAGAGC", 10)),       # HA
      5: (("GGGGGGAGCAGAAGCACAGC", 6),
          ("CCGGGTTATTAGTAGAAACAACAGC", 10)),       # NP
      6: (("GGGGGGAGCAGAAGCAGAGC", 6),
          ("CCGGGTTATTAGTAGTAACAAGAGC", 10)),       # NA
      7: (("GGGGGGAGCAGAAGCACGCACTT", 6),
          ("CCGGGTTATTAGTAGAAACAACGCACTT", 10)),    # M
      8: (("GGGGGGAGCAGAAGCAGAGGATT", 6),
          ("CCGGGTTATTAGTAGTAACAAGAGGATT", 10)),    # NS
  },
}
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


def best_site(primer, seq, skip, first=False):
  """The position of seq most similar to the primer, without gaps (the
  primer anneals with mismatches), outside the first and last <skip>
  bases (the termini): returns (mismatches, start). The best position
  must be unique, unless <first>: the first one is then returned."""
  n = len(primer)
  scores = []
  for i in range(skip, len(seq) - n - skip + 1):
    d = sum(not matches(p, b) for p, b in zip(primer, seq[i:i + n]))
    scores.append((d, i))
  best = min(scores)
  ties = [i for d, i in scores if d == best[0]]
  if len(ties) > 1 and not first:
    raise SystemExit("several equally good sites: %s (see --first-site)"
                     % ties)
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


def pair_of(primers, number):
  """The primers of a segment, as (left, left_labels, right,
  right_labels, left_segment_len, right_segment_len): the left primer as
  it starts a read of the plus strand, the right one as it ends it (its
  reverse complement), each with the origin of each of its bases (t
  tail, p the sequence of the segment)."""
  (left, left_tail), (right, right_tail) = primers[number]
  left_labels = "t" * left_tail + "p" * (len(left) - left_tail)
  # the right primer as it is found at the end of a read: its tail is
  # then at the end, its segment bases at the start
  right_labels = "p" * (len(right) - right_tail) + "t" * right_tail
  return (left, left_labels, revcomp(right), right_labels,
          len(left) - left_tail, len(right) - right_tail)


def amplicon(name, seq, is_ha, site, internal_share, pair, rng,
             wrong_pair=None, wrong_pair_share=0.0):
  """One amplicon of the segment <seq>, as (sequence, labels, kind): full,
  or started at the internal site <site> (a share of those of HA), or
  primed at the segment terminus by the pair of other segments
  (<wrong_pair>, a share), or a primer dimer; possibly with a start or
  an end inside a primer, or broken in the middle."""
  left, left_primer_labels, rc13, rc13_labels, left_len, right_len = pair
  p12 = resolve(left, rng)
  kind = "full"
  r = rng.random()
  # (the share of the amplicons primed by the pair of other segments is
  # drawn only where there is one, so that a simulation without them
  # draws exactly the numbers it drew before)
  wrong = wrong_pair_share > 0 and wrong_pair is not None and \
      rng.random() < wrong_pair_share
  if is_ha and r < internal_share:
    kind = "internal"
    body = seq[site + left_len:len(seq) - right_len]
    left_labels = "i" * len(p12)
  elif wrong:
    # the amplicon began at the terminus of this segment, but was
    # primed by the pair of other segments: the bases of that primer
    # are not the sequence of this segment there
    other_left, _, other_rc13, _, _, _ = wrong_pair
    kind = "wrong_pair"
    p12 = resolve(other_left, rng)
    body = seq[left_len:len(seq) - right_len]
    left_labels = "i" * len(p12)
    rc13, rc13_labels = other_rc13, "i" * len(other_rc13)
  else:
    body = seq[left_len:len(seq) - right_len]
    left_labels = left_primer_labels
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
  return amp, labels, kind


def main():
  parser = argparse.ArgumentParser(
      description="Simulate amplicon reads of influenza with a truth.")
  parser.add_argument("genome")
  parser.add_argument("outdir")
  parser.add_argument(
      "--virus", default="a", choices=sorted(PRIMERS),
      help="the type of influenza the genome is, which chooses the "
           "primers of each segment (default: a)")
  parser.add_argument("--reads-per-segment", type=int, default=400)
  parser.add_argument(
      "--internal-share", type=float, default=0.6,
      help="share of the HA reads started at the internal site "
           "(default 0.6: enough for the artefact to reach a majority "
           "consensus when nothing is trimmed)")
  parser.add_argument("--seed", type=int, default=1)
  parser.add_argument(
      "--error-rate", type=float, default=0.04,
      help="mean of the error rate drawn for each read (normal, standard "
           "deviation half the mean, between 0.005 and 0.15, or 3 times "
           "the mean if higher; default 0.04)")
  parser.add_argument(
      "--barcodes", type=int, default=0, choices=range(0, 5),
      help="native barcoding: each amplicon gets the adapter of one of the "
           "first N barcodes (NB01-NB04) at both ends (default 0: none)")
  parser.add_argument(
      "--segment-reads", default="",
      help="the reads of some segments, as <segment number>:<reads>, "
           "comma-separated (default: --reads-per-segment for all)")
  parser.add_argument(
      "--shuffle", action="store_true",
      help="write the reads in random order (default: by segment)")
  parser.add_argument(
      "--first-site", action="store_true",
      help="if several positions of HA are equally similar to the "
           "primer, use the first as internal priming site (default: "
           "stop, the site must be unique)")
  parser.add_argument(
      "--joined-share", type=float, default=0.0,
      help="share of the reads which are two amplicons joined end to end, "
           "the second one of any segment (default 0: none)")
  parser.add_argument(
      "--wrong-pair-share", type=float, default=0.0,
      help="share of the amplicons primed, at the terminus of their "
           "segment, by the primer pair of other segments (influenza B, "
           "whose primers are not the same in every segment): their "
           "primer bases are not the sequence of the segment there and "
           "must all be removed (default 0: none)")
  parser.add_argument(
      "--read-prefix", default="",
      help="a prefix of the read names (default: none)")
  args = parser.parse_args()
  rng = random.Random(args.seed)
  segment_reads = {}
  for item in filter(None, args.segment_reads.split(",")):
    number, _, reads = item.partition(":")
    segment_reads[int(number)] = int(reads)
  genome = read_fasta(args.genome)
  os.makedirs(args.outdir, exist_ok=True)

  # the segments, named by their number (segment 4 is HA)
  segments = [(name, seq) for name, seq, _ in genome]
  numbers = {}
  for name, _, header in genome:
    m = re.search(r"segment (\d)", header)
    if m:
      numbers[name] = int(m.group(1))
  ha = [(n, s) for n, s, h in genome if "segment 4" in h]
  if len(ha) != 1:
    raise SystemExit("one HA record (\"segment 4\" in its header) is "
                     "needed in the genome")
  ha_name, ha_seq = ha[0]
  primers = PRIMERS[args.virus]
  missing = [n for n in sorted(set(numbers.values())) if n not in primers]
  if missing:
    raise SystemExit("no primers of influenza %s for segment %s"
                     % (args.virus.upper(), missing[0]))
  # the internal priming site: the position of the HA most similar to
  # the primer of its left end, without its tail
  ha_left, ha_left_tail = primers[4][0]
  ha_primer = ha_left[ha_left_tail:]
  dist, site = best_site(ha_primer, ha_seq, 50, args.first_site)
  with open(os.path.join(args.outdir, "site.tsv"), "w") as f:
    f.write("segment\tstart\treplaced\tprimer\tmismatches\n")
    f.write("%s\t%d\t%s\t%s\t%d\n"
            % (ha_name, site, ha_seq[site:site + len(ha_primer)],
               ha_primer, dist))
  with open(os.path.join(args.outdir, "genome.fasta"), "w") as f:
    for _, seq, header in genome:
      f.write(">%s\n%s\n" % (header, seq))

  if args.barcodes > 0:
    with open(os.path.join(args.outdir, "barcodes.fasta"), "w") as f:
      for b in range(args.barcodes):
        f.write(">NB%02d\n%s\n" % (b + 1,
                                    NB_FRONT + NB_BARCODES[b] + NB_REAR))
  # (the primers of influenza A are the same for every segment, so a
  # record whose header does not number its segment can still be
  # simulated; those of influenza B cannot)
  pairs = {}
  for name, _ in segments:
    number = numbers.get(name)
    if number is None:
      if args.virus != "a":
        raise SystemExit("the header of %s does not say which segment it "
                         "is, which the primers of influenza %s need"
                         % (name, args.virus.upper()))
      number = 1
    pairs[name] = pair_of(primers, number)
    numbers.setdefault(name, number)
  # for --wrong-pair-share: for each segment, the pair of a segment
  # whose primers differ from its own (none where every pair is the
  # same, as in influenza A)
  wrong_pairs = {}
  for name, _ in segments:
    other = next((p for n, p in pairs.items()
                  if p[0] != pairs[name][0]), None)
    if other is not None:
      wrong_pairs[name] = other
  records = []
  n = 0
  for name, seq in segments:
    for k in range(segment_reads.get(numbers.get(name),
                                     args.reads_per_segment)):
      amp, labels, kind = amplicon(
          name, seq, name == ha_name, site, args.internal_share,
          pairs[name], rng, wrong_pairs.get(name), args.wrong_pair_share)
      n_amplicons = 1
      if args.joined_share > 0 and rng.random() < args.joined_share:
        # a second amplicon, of any segment, joined to the first one: each
        # in its own orientation, 0 to 20 other bases between them
        other, other_seq = rng.choice(segments)
        amp2, labels2, _ = amplicon(other, other_seq, False, site, 0.0,
                                    pairs[other], rng)
        if rng.random() < 0.5:
          amp, labels = revcomp(amp), labels[::-1]
        if rng.random() < 0.5:
          amp2, labels2 = revcomp(amp2), labels2[::-1]
        gap = rng.randint(0, 20)
        amp += "".join(rng.choice("ACGT") for _ in range(gap)) + amp2
        labels += "f" * gap + labels2
        kind = "joined"
        n_amplicons = 2
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
      rate = min(max(rng.gauss(args.error_rate, args.error_rate / 2),
                     0.005), max(0.15, 3 * args.error_rate))
      read, read_labels = add_errors(amp, labels, rate, rng)
      qual = "".join(chr(33 + rng.randint(8, 30)) for _ in read)
      n += 1
      rname = "%sr%05d" % (args.read_prefix, n)
      records.append(("@%s\n%s\n+\n%s\n" % (rname, read, qual),
                      "%s\t%s\t%s\t%s\t%s\t%s\t%d\n"
                      % (rname, name, kind, strand, read_labels, barcode,
                         n_amplicons)))
  if args.shuffle:
    # (a generator of its own: the reads are the same with and without)
    random.Random(args.seed).shuffle(records)
  fq = gzip.open(os.path.join(args.outdir, "reads.fastq.gz"), "wt")
  truth = gzip.open(os.path.join(args.outdir, "truth.tsv.gz"), "wt")
  truth.write("read\tsegment\tkind\tstrand\tlabels\tbarcode\tamplicons\n")
  for fastq_record, truth_record in records:
    fq.write(fastq_record)
    truth.write(truth_record)
  fq.close()
  truth.close()


if __name__ == "__main__":
  main()
