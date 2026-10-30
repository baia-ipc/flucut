# (c) 2026 Giorgio Gonnella
"""Make reads with the Nanopore adapters, from simulated reads.

For the program tests of the adapter search (testdata/adapters.fastq):
each read of the input (testdata/reads.fastq, see
scripts/simulate-reads.py) gets the adapter of native barcode NB01 at
its start and the reverse complement of the bottom adapter of NB01 at
its end (the sequences of Dorado, see src/adapters.h); in addition, the
outer bases of its primers are removed, as by a break of the molecule:

  <name>_broken      the first 10 bases removed (the read starts inside
                     the primer, or after it), then the adapters added
  <name>_lsk         the ligation adapter (without barcode) at the
                     start, the bottom adapter of NB01 at the end, no
                     base removed

The qualities of the added bases are "I". Only the first <n> reads are
used (default 4).

Usage: python3 -I scripts/add-adapters.py <reads.fastq> [<n>] > out.fastq
"""

import sys

NB01 = "CACAAAGACACCGACAACTTTCTT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
  return seq.translate(COMPLEMENT)[::-1]


# the adapters (src/adapters.h)
NB_START = "ATTGCTAAGGTTAA" + NB01 + "CAGCACCT"
NB_END = revcomp("ATTGCTAAGGTTAA" + NB01 + "CAGCACC")
LSK_START = "CCTGTACTTCGTTCAGTTACGTATTGC"


def main():
  if len(sys.argv) not in (2, 3):
    sys.exit(__doc__)
  n = int(sys.argv[2]) if len(sys.argv) == 3 else 4
  lines = open(sys.argv[1]).read().split("\n")
  for i in range(n):
    title, seq, _, qual = lines[4 * i:4 * i + 4]
    name = title[1:].split()[0]
    for suffix, start, body, bqual in (
        ("broken", NB_START, seq[10:], qual[10:]),
        ("lsk", LSK_START, seq, qual)):
      read = start + body + NB_END
      q = "I" * len(start) + bqual + "I" * len(NB_END)
      print("@%s_%s simulated, with adapters\n%s\n+\n%s"
            % (name, suffix, read, q))


if __name__ == "__main__":
  main()
