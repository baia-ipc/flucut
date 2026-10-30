# (c) 2026 Giorgio Gonnella
"""Count the codons of a site in the reads, by where the reads start.

The reads of a SAM stream (standard input, e.g. samtools view of a
region) that cover the codon at <codon> (1-based, on the reference
<ref>) are classified by the start of their alignment, relative to the
primer-like site starting at <site> (1-based):

  upstream   the alignment starts at least 10 bases before the site: the
             read comes from an amplicon which started before the site
             (e.g. at the segment terminus): its bases at the codon are
             the template
  at_site    the alignment starts within 3 bases of the site: the read
             can come from an amplicon started at the site by internal
             priming: its bases at the codon can be those of the primer
  other      any other start (e.g. a few bases before the site, or
             soft-clipped bases at the start)

Usage: samtools view <bam> <ref>:<codon>-<codon+2> | \\
    python3 -I benchmark/site_reads.py <ref> <site> <codon> <label>
Output: TSV lines label, class, codon, reads (the classes upstream,
at_site, other; the codons by decreasing number of reads).
"""

import collections
import re
import sys


def read_codon(pos, cigar, seq, codon):
  """The read bases at codon..codon+2 (None if not all aligned)."""
  r, q, got = pos, 0, {}
  for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
    n = int(n)
    if op in "M=X":
      for k in range(n):
        if codon <= r + k <= codon + 2:
          got[r + k] = seq[q + k]
      r += n
      q += n
    elif op in "IS":
      q += n
    elif op in "DN":
      r += n
  if len(got) < 3:
    return None
  return "".join(got[codon + k] for k in range(3))


def main():
  if len(sys.argv) != 5:
    sys.exit(__doc__)
  ref, site, codon, label = sys.argv[1], int(sys.argv[2]), \
      int(sys.argv[3]), sys.argv[4]
  counts = collections.Counter()
  for line in sys.stdin:
    f = line.split("\t")
    if f[2] != ref or int(f[1]) & 0x904:
      continue
    pos, cigar, seq = int(f[3]), f[5], f[9]
    c = read_codon(pos, cigar, seq, codon)
    if c is None:
      continue
    if pos <= site - 10:
      cls = "upstream"
    elif abs(pos - site) <= 3:
      cls = "at_site"
    else:
      cls = "other"
    counts[(cls, c)] += 1
  for cls in ("upstream", "at_site", "other"):
    items = sorted(((n, c) for (k, c), n in counts.items() if k == cls),
                   reverse=True)
    for n, c in items:
      print("%s\t%s\t%s\t%d" % (label, cls, c, n))


if __name__ == "__main__":
  main()
