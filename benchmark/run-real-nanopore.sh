#!/bin/bash
# (c) 2026 Giorgio Gonnella
#
# The I351K artefact in a public Nanopore run of H5N1 whose HA consensus
# in GenBank has AAA (Lys) at the codon of the Uni12-like site (the
# Nanopore counterpart of benchmark/run-real-illumina.sh; the reads are
# whole amplicons, single-end, so the default platform of flucut applies
# and IRMA is run as FLU-minion):
#
# 1. the reads (single file, from ENA) and the 8 segments of the
#    consensus of the isolate (GenBank) are downloaded;
# 2. the reads are trimmed by each tool (timed), the others with the
#    primers of the uni12-mbt13 scheme:
#    - flucut, in three or four configurations: default options (the
#      scheme chosen from the reads, -s auto; the
#      reads classified against the built-in references); -r off (no
#      references); -r with the 8 segments of the isolate (coding
#      sequences in most records: the terminal read ends are then not
#      located); and, if REFERENCES is set, -r with that file (complete
#      segments, with the untranslated regions);
#    - cutadapt, the 5' Uni12 and 3' Uni13 (reverse complement), both
#      orientations (--revcomp);
#    - BBDuk, the primers as k-mers at both ends (ktrim=l, then ktrim=r);
#    - the trimmer of IRMA-core, the primers as k-mers (window 120);
#    - dorado trim, the primers as a primer file (if DORADO is set);
# 3. the codon is counted in the reads of each tool by where they map
#    (site_reads.py, minimap2 -ax map-ont): the reads starting upstream
#    of the site come from amplicons started before it, and show the
#    template;
# 4. IRMA FLU-minion assembles the reads of each tool, and the consensus
#    sequences are compared with that of the untrimmed reads.
#
# Usage: benchmark/run-real-nanopore.sh <run> <accessions> <output dir> \
#          [<download dir>]
# <accessions>: the 8 segments of the isolate, comma-separated. The reads
# and the segments are downloaded to <download dir> (default: <output
# dir>/download; nothing is downloaded again if it is there), which
# holds nothing else: no program is run from inside it. The programs are
# given in the environment: FLUCUT, IRMA, SAMTOOLS, MINIMAP2 (required);
# IRMA_CORE, BBDUK, DORADO, REFERENCES (optional); cutadapt in the PATH
# (optional). IRMA_JOBS: the IRMA runs at a time (default 1);
# MEASURE_CPUS: see measure.py. Output: <output dir>/speed.tsv (time and
# memory), bases.tsv (reads and bases written), codons.tsv,
# consensus.tsv.

set -euo pipefail
run=$1
accessions=$2
out=$3
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
measure="python3 -I $here/measure.py"
uni12=AGCRAAAGCAGG
uni12a=AGCAAAAGCAGG
uni12g=AGCGAAAGCAGG
mbt13=ACGCGTGATCAGTAGAAACAAGG
rc_mbt13=CCTTGTTTCTACTGATCACGCGT
mkdir -p "$out"
out=$(cd "$out" && pwd)
dl=${4:-$out/download}
mkdir -p "$dl"
dl=$(cd "$dl" && pwd)
isolate=$dl/isolate.fasta
cd "$out"

# 1. the data
if [ ! -s "$dl/$run.fastq.gz" ]; then
  url=$(curl -sS --fail --retry 3 "https://www.ebi.ac.uk/ena/portal/api/\
filereport?accession=$run&result=read_run&fields=fastq_ftp" \
    | tail -1 | cut -f2 || true)
  case $url in
    *ftp.sra.ebi.ac.uk*) ;;
    *) d=ftp.sra.ebi.ac.uk/vol1/fastq/${run:0:6}/0${run: -2}/$run
       url=$d/${run}_1.fastq.gz ;;
  esac
  curl -sS --fail --retry 5 --retry-all-errors --retry-delay 10 \
    -o "$dl/$run.fastq.gz" "https://${url%%;*}"
fi
if [ ! -s "$isolate" ]; then
  curl -sS --fail --retry 5 --retry-all-errors --retry-delay 10 \
    -o "$isolate" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/\
efetch.fcgi?db=nuccore&id=$accessions&rettype=fasta&retmode=text"
fi
[ -s "$dl/download_date.txt" ] || date -u +%Y-%m-%d > "$dl/download_date.txt"
reads=$dl/$run.fastq.gz
# the HA (segment 4) and its Uni12-like site (12 bases, 1-based start;
# the codon is its bases 5-7), searched near position 1050
read -r ha site < <(python3 -I - "$isolate" <<'EOF'
import re, sys
seqs, name = {}, None
for line in open(sys.argv[1]):
  line = line.strip()
  if line.startswith(">"):
    name = line[1:].split()[0]
    seqs[name] = [line]
  else:
    seqs[name].append(line.upper())
for name, lines in seqs.items():
  if "segment 4" in lines[0]:
    ha = "".join(lines[1:])
    m = re.search("[AG]GC[AGT]A[AT]AGCAGG", ha[1000:1100])
    print(name, m.start() + 1001)
EOF
)
codon=$(( site + 4 ))
printf 'run\t%s\nHA\t%s\nsite\t%s\ncodon\t%s\n' $run $ha $site $codon \
  > site.tsv

# 2. the reads trimmed by each tool (<tool>.fastq.gz), timed
printf 'run\ttool\twall_s\tcpu_s\tpeak_mb\n' > speed.tsv
timed() { $measure "$1" "$2" | sed "s/^/$run\t/" >> speed.tsv; }
tools="untrimmed flucut flucut-r-off flucut-isolate"
ln -sf $reads untrimmed.fastq.gz
rm -rf flucut flucut-r-off flucut-isolate flucut-references
timed flucut "$FLUCUT -q -o flucut $reads"
timed flucut-r-off "$FLUCUT -q -r off -o flucut-r-off $reads"
timed flucut-isolate "$FLUCUT -q -r $isolate -o flucut-isolate $reads"
if [ -n "${REFERENCES:-}" ]; then
  tools="$tools flucut-references"
  timed flucut-references "$FLUCUT -q -r $REFERENCES -o flucut-references \
    $reads"
fi
for t in $tools; do
  [ $t = untrimmed ] || ln -sf $t/$run.fastq.gz $t.fastq.gz
done
if command -v cutadapt > /dev/null; then
  tools="$tools cutadapt"
  timed cutadapt "cutadapt -j 4 -e 0.2 --revcomp -n 2 -m 1 -g $uni12 \
    -a $rc_mbt13 -o cutadapt.fastq.gz $reads > cutadapt.log"
fi
if [ -n "${BBDUK:-}" ]; then
  tools="$tools bbduk"
  printf '>uni12A\n%s\n>uni12G\n%s\n>mbtuni13\n%s\n' $uni12a $uni12g \
    $mbt13 > primers.fasta
  common="ref=primers.fasta k=12 mink=8 hdist=2 rcomp=t threads=4 -Xmx2g ow=t"
  timed bbduk "$BBDUK in=$reads out=bbduk-l.fastq.gz ktrim=l \
    restrictleft=100 $common 2> bbduk-l.log && $BBDUK in=bbduk-l.fastq.gz \
    out=bbduk.fastq.gz ktrim=r restrictright=100 minlen=1 $common \
    2> bbduk-r.log"
fi
if [ -n "${IRMA_CORE:-}" ]; then
  tools="$tools irma-core"
  [ -s primers.fasta ] || printf '>uni12A\n%s\n>uni12G\n%s\n>mbtuni13\n%s\n' \
    $uni12a $uni12g $mbt13 > primers.fasta
  timed irma-core "$IRMA_CORE trimmer $reads -o irma-core.fastq \
    -P primers.fasta --p-kmer-length 12 --p-fuzzy --p-restrict 120 \
    2> irma-core.log && gzip -f irma-core.fastq"
fi
if [ -n "${DORADO:-}" ]; then
  tools="$tools dorado"
  dp=dorado-primers.fasta
  : > "$dp"
  for form in A:$uni12a G:$uni12g; do
    printf '>uni12%s_front\tet:Z:primer\tsk:Z:ANY\n%s\n' "${form%%:*}" \
      "${form#*:}" >> "$dp"
    printf '>uni12%s_rear\tet:Z:primer\tsk:Z:ANY\n%s\n' "${form%%:*}" \
      "$mbt13" >> "$dp"
  done
  timed dorado "$DORADO trim -k SQK-LSK114 --primer-sequences $dp \
    --emit-fastq -t 4 $reads 2> dorado.log | gzip > dorado.fastq.gz"
fi

# 3. the codon in the reads of each tool, and the bases written
printf 'run\treads\tclass\tcodon\tcount\n' > codons.tsv
printf 'run\ttool\treads\tbases\n' > bases.tsv
for t in $tools; do
  printf '%s\t%s\t%s\n' $run $t "$(gzip -dc ${t}.fastq.gz | awk \
    'NR % 4 == 2 { n++; b += length($0) } END { print n "\t" b }')" \
    >> bases.tsv
  "$MINIMAP2" -ax map-ont -t 4 "$isolate" ${t}.fastq.gz 2> /dev/null \
    | "$SAMTOOLS" sort -o $t.bam -
  "$SAMTOOLS" index $t.bam
  "$SAMTOOLS" view $t.bam $ha:$codon-$((codon + 2)) \
    | python3 -I "$here/site_reads.py" $ha $site $codon $t \
    | sed "s/^/$run\t/" >> codons.tsv
done

# 4. the consensus of IRMA (FLU-minion)
for t in $tools; do
  rm -rf irma_$t
  "$IRMA" FLU-minion ${t}.fastq.gz irma_$t > irma_$t.log 2>&1 || true &
  while [ "$(jobs -rp | wc -l)" -ge "${IRMA_JOBS:-1}" ]; do wait -n; done
done
wait
printf 'run\treads\tsegment\tlength\tcodon\tdifferences_to_untrimmed\n' \
  > consensus.tsv
python3 -I - $run $site $tools >> consensus.tsv <<'EOF'
import glob, os, sys
import edlib
run, site, tools = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
def read(path):
  return "".join(l.strip() for l in open(path) if not l.startswith(">"))
base = sorted(glob.glob("irma_untrimmed/A_*.fasta"))
for f in base:
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
      i = s.find(ref[site - 11:site - 1])
      codon = s[i + 14:i + 17] if i >= 0 else "not found"
    print("%s\t%s\t%s\t%d\t%s\t%d" % (run, t, seg, len(s), codon, d))
EOF
cat site.tsv speed.tsv bases.tsv codons.tsv consensus.tsv
