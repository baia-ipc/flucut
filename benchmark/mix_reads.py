# (c) 2026 Giorgio Gonnella
"""Put the reads of several simulations together, in random order.

For samples with several strains (a mixed infection): each strain is
simulated by itself with benchmark/simulate.py, with the number of
reads of its share of the sample and its own --read-prefix (the read
names must differ), and the simulations are joined here.

Writes reads.fastq.gz and truth.tsv.gz to <outdir>, as simulate.py does
(the truth of each read is that of its simulation).

Usage:
  python3 -I benchmark/mix_reads.py <outdir> <seed> <simulation dir>...
"""

import gzip
import os
import random
import sys


def main():
  if len(sys.argv) < 4:
    sys.exit(__doc__)
  outdir, seed, simdirs = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
  records, header, names = [], "", set()
  for simdir in simdirs:
    with gzip.open(simdir + "/reads.fastq.gz", "rt") as fq, \
         gzip.open(simdir + "/truth.tsv.gz", "rt") as truth:
      header = truth.readline()
      for line in truth:
        read = "".join(fq.readline() for _ in range(4))
        name = line.split("\t")[0]
        if name in names:
          sys.exit("%s: read name in several simulations (use "
                   "--read-prefix)" % name)
        names.add(name)
        records.append((read, line))
  random.Random(seed).shuffle(records)
  os.makedirs(outdir, exist_ok=True)
  with gzip.open(outdir + "/reads.fastq.gz", "wt") as fq, \
       gzip.open(outdir + "/truth.tsv.gz", "wt") as truth:
    truth.write(header)
    for read, line in records:
      fq.write(read)
      truth.write(line)


if __name__ == "__main__":
  main()
