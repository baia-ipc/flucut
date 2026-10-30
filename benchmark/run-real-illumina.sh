#!/bin/bash
# (c) 2026 Giorgio Gonnella
#
# The I348K artefact in a public Illumina run of H5N1 whose HA consensus
# in GenBank has AAA (Lys) at the codon of the Uni12-like site (the runs
# were found by a scan of the public H5N1 HA records with that codon, see
# benchmark/RESULTS.md):
#
# 1. the reads (paired, from ENA) and the 8 segments of the consensus of
#    the isolate (GenBank) are downloaded;
# 2. the reads are trimmed by each tool, R1 and R2 as two inputs:
#    - flucut -p illumina, with the 8 segments as references;
#    - cutadapt, the 5' primers (Uni12, Uni13) not anchored, and anchored
#      (^; the reads start inside the tail of the primers, which is not
#      given: up to 2 errors);
#    - BBDuk, the primers as k-mers (k 12, hdist 1) in the first 30
#      bases (ktrim=l);
#    - the trimmer of IRMA-core, the primers as k-mers in the first 30
#      bases;
# 3. the codon is counted in the reads of each tool by where they start
#    (site_reads.py): the reads starting upstream of the site come from
#    amplicons started before it, and show the template;
# 4. IRMA FLU assembles the reads of each tool, and the consensus
#    sequences are compared with that of the untrimmed reads (the codon,
#    and the differences in the 8 segments).
#
# Usage: benchmark/run-real-illumina.sh <run> <accessions> <output dir>
# <accessions>: the 8 segments of the isolate, comma-separated. The
# programs are given in the environment: FLUCUT, IRMA, SAMTOOLS, MINIMAP2
# (required); IRMA_CORE, BBDUK (optional); cutadapt in the PATH
# (optional). Output: <output dir>/speed.tsv (time and memory of the
# trimming, both files), bases.tsv (reads and bases written),
# codons.tsv, consensus.tsv.

set -euo pipefail
run=$1
accessions=$2
out=$3
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
uni12=AGCRAAAGCAGG
uni12a=AGCAAAAGCAGG
uni12g=AGCGAAAGCAGG
uni13=AGTAGAAACAAGG
mkdir -p "$out"
cd "$out"

# 1. the data
if [ ! -s ${run}_2.fastq.gz ]; then
  urls=$(curl -sS --fail "https://www.ebi.ac.uk/ena/portal/api/filereport?\
accession=$run&result=read_run&fields=fastq_ftp" | tail -1 | cut -f2)
  for url in ${urls//;/ }; do
    curl -sS --fail -o "$(basename "$url")" "https://$url"
  done
fi
if [ ! -s isolate.fasta ]; then
  curl -sS --fail -o isolate.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/\
eutils/efetch.fcgi?db=nuccore&id=$accessions&rettype=fasta&retmode=text"
fi
[ -s download_date.txt ] || date -u +%Y-%m-%d > download_date.txt
# the HA (segment 4) and its Uni12-like site (12 bases, 1-based start;
# the codon is its bases 5-7), searched near position 1050
read -r ha site < <(python3 -I - <<'EOF'
import re
seqs, name = {}, None
for line in open("isolate.fasta"):
  line = line.strip()
  if line.startswith(">"):
    name = line[1:].split()[0]
    seqs[name] = [line]
  else:
    seqs[name].append(line.upper())
for name, lines in seqs.items():
  if "segment 4" in lines[0]:
    ha = "".join(lines[1:])
    m = re.search("GGC[GT][AT]A[AT]GCAGG", ha[1000:1100])
    print(name, m.start() + 1001)
EOF
)
codon=$(( site + 4 ))
printf 'run\t%s\nHA\t%s\nsite\t%s\ncodon\t%s\n' $run $ha $site $codon \
  > site.tsv

# 2. the reads trimmed by each tool (<tool>_1.fastq.gz, <tool>_2...),
#    timed (measure.py: wall s, CPU s, peak MB; both files) in speed.tsv
measure="python3 -I $here/measure.py"
printf 'run\ttool\twall_s\tcpu_s\tpeak_mb\n' > speed.tsv
timed() { $measure "$1" "$2" | sed "s/^/$run\t/" >> speed.tsv; }
tools="untrimmed flucut"
for i in 1 2; do
  ln -sf ${run}_$i.fastq.gz untrimmed_$i.fastq.gz
done
rm -rf flucut
timed flucut "$FLUCUT -q -p illumina -r isolate.fasta -o flucut \
  ${run}_1.fastq.gz ${run}_2.fastq.gz"
for i in 1 2; do
  ln -sf flucut/${run}_$i.fastq.gz flucut_$i.fastq.gz
done
if command -v cutadapt > /dev/null; then
  tools="$tools cutadapt cutadapt-anchored"
  timed cutadapt "for i in 1 2; do cutadapt -j 4 -e 0.2 -m 1 \
    -g $uni12 -g $uni13 -o cutadapt_\$i.fastq.gz ${run}_\$i.fastq.gz \
    > cutadapt_\$i.log; done"
  timed cutadapt-anchored "for i in 1 2; do cutadapt -j 4 -e 0.2 -m 1 \
    -g ^$uni12a -g ^$uni12g -g ^$uni13 -o cutadapt-anchored_\$i.fastq.gz \
    ${run}_\$i.fastq.gz > cutadapt-anchored_\$i.log; done"
fi
if [ -n "${BBDUK:-}" ]; then
  tools="$tools bbduk"
  printf '>uni12A\n%s\n>uni12G\n%s\n>uni13\n%s\n' $uni12a $uni12g $uni13 \
    > primers.fasta
  timed bbduk "for i in 1 2; do $BBDUK in=${run}_\$i.fastq.gz \
    out=bbduk_\$i.fastq.gz ref=primers.fasta k=12 mink=8 hdist=1 \
    ktrim=l restrictleft=30 minlen=1 threads=4 -Xmx2g ow=t \
    2> bbduk_\$i.log; done"
fi
if [ -n "${IRMA_CORE:-}" ]; then
  tools="$tools irma-core"
  [ -s primers.fasta ] || printf '>uni12A\n%s\n>uni12G\n%s\n>uni13\n%s\n' \
    $uni12a $uni12g $uni13 > primers.fasta
  timed irma-core "for i in 1 2; do $IRMA_CORE trimmer \
    ${run}_\$i.fastq.gz -o irma-core_\$i.fastq -P primers.fasta \
    --p-kmer-length 12 --p-fuzzy --p-restrict 30 2> irma-core_\$i.log; \
    gzip -f irma-core_\$i.fastq; done"
fi

# 3. the codon in the reads of each tool, and the bases written
printf 'run\treads\tclass\tcodon\tcount\n' > codons.tsv
printf 'run\ttool\treads\tbases\n' > bases.tsv
for t in $tools; do
  printf '%s\t%s\t%s\n' $run $t "$(gzip -dc ${t}_1.fastq.gz \
    ${t}_2.fastq.gz | awk 'NR % 4 == 2 { n++; b += length($0) }
    END { print n "\t" b }')" >> bases.tsv
  "$MINIMAP2" -ax sr -t 4 isolate.fasta ${t}_1.fastq.gz ${t}_2.fastq.gz \
    2> /dev/null | "$SAMTOOLS" sort -o $t.bam -
  "$SAMTOOLS" index $t.bam
  "$SAMTOOLS" view $t.bam $ha:$codon-$((codon + 2)) \
    | python3 -I "$here/site_reads.py" $ha $site $codon $t \
    | sed "s/^/$run\t/" >> codons.tsv
done

# 4. the consensus of IRMA
for t in $tools; do
  rm -rf irma_$t
  "$IRMA" FLU ${t}_1.fastq.gz ${t}_2.fastq.gz irma_$t > irma_$t.log 2>&1
done
printf 'run\treads\tsegment\tlength\tcodon\tdifferences_to_untrimmed\n' \
  > consensus.tsv
python3 -I - $run $site $tools >> consensus.tsv <<'EOF'
import glob, os, sys
import edlib
run, site, tools = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
def read(path):
  return "".join(l.strip() for l in open(path) if not l.startswith(">"))
for f in sorted(glob.glob("irma_untrimmed/A_*.fasta")):
  seg = os.path.basename(f)[2:-6]
  ref = read(f).upper()
  for t in tools:
    g = "irma_%s/%s" % (t, os.path.basename(f))
    if not os.path.exists(g):
      print("%s\t%s\t%s\tNA\tNA\tNA" % (run, t, seg))
      continue
    s = read(g).upper()
    d = edlib.align(s, ref)["editDistance"]
    codon = "NA"
    if seg.startswith("HA"):
      # the codon, after the 10 bases before the site (in the isolate
      # consensus and in that of the untrimmed reads, the same length)
      i = s.find(ref[site - 11:site - 1])
      codon = s[i + 14:i + 17] if i >= 0 else "not found"
    print("%s\t%s\t%s\t%d\t%s\t%d" % (run, t, seg, len(s), codon, d))
EOF
cat site.tsv speed.tsv bases.tsv codons.tsv consensus.tsv
