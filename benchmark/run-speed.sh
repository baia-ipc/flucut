#!/bin/bash
# (c) 2026 Giorgio Gonnella
#
# Time and memory of the trimming tools (benchmark/measure.py), with the
# commands of benchmark/run-tools.sh (see there), on the reads of a
# simulation repeated <copies> times (the read names made unique), with 4
# threads where the tool has an option for it (and flucut also with 1).
# The reads written are not scored here (see run-tools.sh). Barbell needs
# barcoded reads: it is run on the reads of <barcoded simulation>, if
# given, repeated as well (the other tools are run on them too, then).
#
# Usage: benchmark/run-speed.sh <simulation dir> <output dir> [<copies>]
# The programs are given in the environment, as for run-tools.sh.
# Output: <output dir>/speed.tsv (tool, wall s, CPU s, peak MB).

set -euo pipefail
sim=$1
out=$2
copies=${3:-10}
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
measure="python3 -I $here/measure.py"
uni12=AGCRAAAGCAGG
uni12a=AGCAAAAGCAGG
uni12g=AGCGAAAGCAGG
mbt13=ACGCGTGATCAGTAGAAACAAGG
rc_mbt13=CCTTGTTTCTACTGATCACGCGT
mkdir -p "$out"
reads=$out/reads.fastq.gz
for (( i = 0; i < copies; i++ )); do
  gzip -dc "$sim/reads.fastq.gz" |
    awk -v c=$i 'NR % 4 == 1 { sub(/^@/, "@c" c "_") } 1'
done | gzip > "$reads"
primers=$out/primers.fasta
printf '>uni12A\n%s\n>uni12G\n%s\n>mbtuni13\n%s\n' "$uni12a" "$uni12g" \
  "$mbt13" > "$primers"
table=$out/speed.tsv
printf 'tool\twall_s\tcpu_s\tpeak_mb\n' > "$table"
run() { $measure "$1" "$2" >> "$table"; }

run flucut-j1 "$FLUCUT -q -j 1 -o $out/flucut1 $reads"
run flucut-j4 "$FLUCUT -q -j 4 -o $out/flucut4 $reads"
run flucut-j4-references "$FLUCUT -q -j 4 -r $sim/genome.fasta \
  -o $out/flucut4r $reads"
if command -v cutadapt > /dev/null; then
  run cutadapt-j4 "cutadapt -j 4 -e 0.2 --revcomp -n 2 -m 1 -g $uni12 \
    -a $rc_mbt13 -o $out/cutadapt.fastq.gz $reads > /dev/null"
fi
if [ -n "${IRMA_CORE:-}" ]; then
  run irma-core-w70 "$IRMA_CORE trimmer $reads -o $out/irma-core.fastq \
    -P $primers --p-kmer-length 12 --p-fuzzy --p-restrict 70 -n 1 \
    2> /dev/null && gzip -f $out/irma-core.fastq"
fi
if [ -n "${BBDUK:-}" ]; then
  common="ref=$primers k=12 mink=8 hdist=2 rcomp=t threads=4 -Xmx2g ow=t"
  run bbduk-hdist2-t4 "$BBDUK in=$reads out=$out/bbduk-l.fastq.gz ktrim=l \
    restrictleft=100 $common 2> /dev/null && $BBDUK \
    in=$out/bbduk-l.fastq.gz out=$out/bbduk.fastq.gz ktrim=r \
    restrictright=100 minlen=1 $common 2> /dev/null"
fi
if [ -n "${PORECHOP:-}" ]; then
  adapters=$out/porechop-adapters.txt
  printf 'uni12A\n%s\n%s\nuni12G\n%s\n%s\n' "$uni12a" "$rc_mbt13" \
    "$uni12g" "$rc_mbt13" > "$adapters"
  run porechop-abi-tuned-t4 "$PORECHOP -i $reads -o $out/porechop.fastq.gz \
    -cap $adapters -ddb -t 4 --extra_end_trim 0 --end_threshold 60 \
    > /dev/null 2>&1"
fi
if [ -n "${DORADO:-}" ]; then
  dp=$out/dorado-primers.fasta
  : > "$dp"
  for form in A:$uni12a G:$uni12g; do
    printf '>uni12%s_front\tet:Z:primer\tsk:Z:ANY\n%s\n' "${form%%:*}" \
      "${form#*:}" >> "$dp"
    printf '>uni12%s_rear\tet:Z:primer\tsk:Z:ANY\n%s\n' "${form%%:*}" \
      "$mbt13" >> "$dp"
  done
  run dorado-trim-t4 "$DORADO trim -k SQK-LSK114 --primer-sequences $dp \
    --emit-fastq -t 4 $reads 2> /dev/null | gzip > $out/dorado.fastq.gz"
fi
if [ -n "${SAMTOOLS:-}" ] && [ -n "${MINIMAP2:-}" ]; then
  genome=$out/genome.fasta
  cp "$sim/genome.fasta" "$genome"
  "$SAMTOOLS" faidx "$genome"
  awk -v OFS='\t' '{ print $1, 0, 12, "uni12_" $1, 0, "+"
                     print $1, $2 - 13, $2, "uni13_" $1, 0, "-" }' \
    "$genome.fai" > "$out/primers.bed"
  # (the mapping, the sorting and the clipping; not the conversion back
  # to FASTQ of run-tools.sh)
  run minimap2-ampliconclip-t4 "$MINIMAP2 -ax map-ont -t 4 $genome $reads \
    2> /dev/null | $SAMTOOLS view -b -F 0x904 - | $SAMTOOLS sort \
    -o $out/mapped.bam - && $SAMTOOLS index $out/mapped.bam && \
    $SAMTOOLS ampliconclip --hard-clip --both-ends -b $out/primers.bed \
    -o $out/clipped.bam $out/mapped.bam 2> /dev/null"
fi
if [ -n "${BARBELL:-}" ] && [ -e "$sim/barcodes.fasta" ]; then
  bb=$out/barbell
  mkdir -p "$bb"
  awk -v p="$uni12" '/^>/ { print; next } { print $0 p }' \
    "$sim/barcodes.fasta" > "$bb/left.fasta"
  awk -v p="$mbt13" '/^>/ { print; next } { print $0 p }' \
    "$sim/barcodes.fasta" > "$bb/right.fasta"
  {
    echo 'Ftag[fw, *, @left(0..250), >>]__Rtag[<<, rc, *, @right(0..250)]'
    echo 'Rtag[fw, *, @left(0..250), >>]__Ftag[<<, rc, *, @right(0..250)]'
    echo 'Ftag[fw, *, @left(0..250), >>]__Ftag[<<, rc, *, @right(0..250)]'
    echo 'Rtag[fw, *, @left(0..250), >>]__Rtag[<<, rc, *, @right(0..250)]'
    echo 'Ftag[fw, *, @left(0..250), >>]'
    echo 'Rtag[fw, *, @left(0..250), >>]'
    echo 'Ftag[<<, rc, *, @right(0..250)]'
    echo 'Rtag[<<, rc, *, @right(0..250)]'
  } > "$bb/filters.txt"
  run barbell-t4 "$BARBELL annotate -q $bb/left.fasta $bb/right.fasta \
    -b Ftag Rtag -i $reads -o $bb/anno.tsv -t 4 > /dev/null 2>&1 && \
    $BARBELL filter -i $bb/anno.tsv -f $bb/filters.txt \
    -o $bb/filtered.tsv > /dev/null 2>&1 && $BARBELL trim \
    -i $bb/filtered.tsv -r $reads -o $bb/trimmed --no-orientation \
    --no-flanks --sort-labels > /dev/null 2>&1"
fi
cat "$table"
