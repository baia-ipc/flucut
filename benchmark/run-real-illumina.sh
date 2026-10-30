#!/bin/bash
# (c) 2026 Giorgio Gonnella
#
# The I348K artefact in a public Illumina run (SRA SRR34913296, H5N1,
# A/cattle/CA/25-020784-009-original/2025, NextSeq 2000, paired; its HA
# consensus in GenBank, PX279898.1, has AAA (Lys) at the codon of the
# Uni12-like site; the run was found by a scan of the public H5N1 HA
# records with that codon, see benchmark/RESULTS.md):
#
# 1. the reads (ENA) and the 8 segments of the isolate (GenBank
#    PX279895.1-PX279902.1) are downloaded;
# 2. the codon is counted in the reads by where they start
#    (site_reads.py): the reads starting upstream of the site come from
#    amplicons started before it, and show the template;
# 3. flucut -p illumina trims the reads (R1 and R2 as two inputs);
# 4. IRMA FLU assembles the reads untrimmed and trimmed, and the
#    consensus sequences are compared (the codon, and all differences).
#
# Usage: benchmark/run-real-illumina.sh <output dir>
# The programs are given in the environment: FLUCUT, IRMA, SAMTOOLS,
# MINIMAP2 (see run-tools.sh). Output: <output dir>/codons.tsv,
# <output dir>/consensus.tsv, and the runs.

set -euo pipefail
out=$1
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
run=SRR34913296
ha=PX279898.1
mkdir -p "$out"
cd "$out"

# 1. the data
for i in 1 2; do
  [ -s ${run}_$i.fastq.gz ] || curl -sS --fail -o ${run}_$i.fastq.gz \
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR349/096/$run/${run}_$i.fastq.gz"
done
if [ ! -s isolate.fasta ]; then
  ids=PX279895.1,PX279896.1,PX279897.1,PX279898.1,PX279899.1
  ids=$ids,PX279900.1,PX279901.1,PX279902.1
  curl -sS --fail -o isolate.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/\
eutils/efetch.fcgi?db=nuccore&id=$ids&rettype=fasta&retmode=text"
fi
date -u +%Y-%m-%d > download_date.txt
# the Uni12-like site of the HA (12 bases, 1-based start; the codon is
# its bases 5-7): searched near position 1050
site=$(python3 -I - "$ha" <<'EOF'
import re, sys
seqs, name = {}, None
for line in open("isolate.fasta"):
  line = line.strip()
  if line.startswith(">"):
    name = line[1:].split()[0]
    seqs[name] = []
  else:
    seqs[name].append(line.upper())
ha = "".join(seqs[sys.argv[1]])
m = re.search("GGC[GT][AT]A[AT]GCAGG", ha[1000:1100])
print(m.start() + 1001)
EOF
)
codon=$(( site + 4 ))
echo "site $site, codon $codon" > site.txt

# 2. the codon in the reads, untrimmed and trimmed (3.)
map() {
  "$MINIMAP2" -ax sr -t 4 isolate.fasta "$2" "$3" 2> /dev/null \
    | "$SAMTOOLS" sort -o "$1.bam" -
  "$SAMTOOLS" index "$1.bam"
}
rm -rf flucut
"$FLUCUT" -q -p illumina -r isolate.fasta -o flucut ${run}_1.fastq.gz \
  ${run}_2.fastq.gz
map untrimmed ${run}_1.fastq.gz ${run}_2.fastq.gz
map trimmed flucut/${run}_1.fastq.gz flucut/${run}_2.fastq.gz
printf 'reads\tclass\tcodon\tcount\n' > codons.tsv
for b in untrimmed trimmed; do
  "$SAMTOOLS" view $b.bam $ha:$codon-$((codon + 2)) \
    | python3 -I "$here/site_reads.py" $ha $site $codon $b >> codons.tsv
done

# 4. the consensus of IRMA
printf 'reads\tsegment\tlength\tcodon\tdifferences_to_untrimmed\n' \
  > consensus.tsv
for b in untrimmed trimmed; do
  rm -rf irma_$b
  if [ $b = untrimmed ]; then
    "$IRMA" FLU ${run}_1.fastq.gz ${run}_2.fastq.gz irma_$b > irma_$b.log 2>&1
  else
    "$IRMA" FLU flucut/${run}_1.fastq.gz flucut/${run}_2.fastq.gz irma_$b \
      > irma_$b.log 2>&1
  fi
done
python3 -I - "$site" >> consensus.tsv <<'EOF'
import glob, os, sys
import edlib
site = int(sys.argv[1])
def read(path):
  return "".join(l.strip() for l in open(path) if not l.startswith(">"))
for f in sorted(glob.glob("irma_untrimmed/A_*.fasta")):
  seg = os.path.basename(f)[2:-6]
  ref = read(f).upper()
  for b in ("untrimmed", "trimmed"):
    g = "irma_%s/%s" % (b, os.path.basename(f))
    if not os.path.exists(g):
      print("%s\t%s\tNA\tNA\tNA" % (b, seg))
      continue
    s = read(g).upper()
    d = edlib.align(s, ref)["editDistance"]
    codon = "NA"
    if seg.startswith("HA"):
      # the codon, after the 10 bases before the site (as in the
      # untrimmed consensus)
      before = ref[site - 11:site - 1]
      i = s.find(before)
      codon = s[i + 14:i + 17] if i >= 0 else "not found"
    print("%s\t%s\t%d\t%s\t%d" % (b, seg, len(s), codon, d))
EOF
cat site.txt codons.tsv consensus.tsv
