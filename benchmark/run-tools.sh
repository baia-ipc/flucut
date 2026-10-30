#!/bin/bash
#
# The comparison of the tools on simulated reads (a directory written by
# simulate.py): each tool trims the reads, the trimmed reads are scored
# against the truth (score_reads.py), and IRMA FLU-minion assembles them,
# as it is and with the reference elongation (SKIP_E=0); the consensus is
# compared with the truth (compare_consensus.py), and the codon of the
# I348K site is read in the plain consensus of HA.
#
# The tools, with the configurations of benchmark/RESULTS.md (each tool
# runs if its program is given, see below):
# - flucut, default options;
# - cutadapt: the two primers non-anchored (-g, -a), --revcomp, -e 0.2;
# - IRMA-core trimmer: the primers as k-mers (k 12, the length of Uni12;
#   one mismatch), windows of 40 (the default of MIRA-NF for custom
#   primers), 70 and 120 bases;
# - BBDuk: two passes (ktrim=l, then ktrim=r), k 12, mink 8, hdist 1 or
#   2, both strands, windows of 100 bases;
# - Porechop_ABI: the primers as custom adapters (-cap, -ddb), with its
#   defaults, and with --extra_end_trim 0 --end_threshold 60 ("tuned");
# - dorado trim: the primers as custom primers (sk:Z:ANY), kit
#   SQK-LSK114; for barcoded reads, dorado demux (SQK-NBD114-24) first;
# - samtools ampliconclip: the reads mapped by minimap2 (map-ont) to the
#   true genome, the primer intervals (the first 12 and the last 13
#   bases of each segment) hard-clipped (--both-ends), the soft-clipped
#   ends removed; unmapped reads are not written;
# - Barbell (barcoded reads only: barcodes.fasta of the simulation):
#   queries made of the barcoding adapter and the primer, the patterns of
#   one or both ends, cut after the left query and before the right one.
#
# Usage: benchmark/run-tools.sh <simdir> <outdir>
# The programs are given in the environment: FLUCUT, IRMA (required);
# IRMA_CORE, BBDUK, PORECHOP, DORADO, SAMTOOLS and MINIMAP2 (together),
# BARBELL (optional); cutadapt is run if it is in the PATH.
# Output: <outdir>/{reads,consensus,site}.tsv, also shown.
#

set -euo pipefail
sim=$1
out=$2
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
reads=$sim/reads.fastq.gz
uni12=AGCRAAAGCAGG
uni12a=AGCAAAAGCAGG
uni12g=AGCGAAAGCAGG
mbt13=ACGCGTGATCAGTAGAAACAAGG
rc_mbt13=CCTTGTTTCTACTGATCACGCGT
barcoded=0
[ -e "$sim/barcodes.fasta" ] && barcoded=1
mkdir -p "$out"
printf 'SKIP_E=0\n' > "$out/elong.sh"
primers=$out/primers.fasta
printf '>uni12A\n%s\n>uni12G\n%s\n>mbtuni13\n%s\n' "$uni12a" "$uni12g" \
  "$mbt13" > "$primers"

# the trimmed reads of each tool, in $out/<tool>.fastq.gz
tools="untrimmed flucut"
cp "$reads" "$out/untrimmed.fastq.gz"
"$FLUCUT" -q -o "$out/flucut" "$reads"
mv "$out/flucut/reads.fastq.gz" "$out/flucut.fastq.gz"

if command -v cutadapt > /dev/null; then
  tools="$tools cutadapt"
  cutadapt -j 1 -e 0.2 --revcomp -n 2 -m 1 -g "$uni12" -a "$rc_mbt13" \
    -o "$out/cutadapt.fastq.gz" "$reads" > "$out/cutadapt.log"
fi

if [ -n "${IRMA_CORE:-}" ]; then
  for w in 40 70 120; do
    tools="$tools irma-core-w$w"
    "$IRMA_CORE" trimmer "$reads" -o "$out/irma-core-w$w.fastq" \
      -P "$primers" --p-kmer-length 12 --p-fuzzy --p-restrict "$w" -n 1 \
      -v 2> "$out/irma-core-w$w.log"
    gzip -f "$out/irma-core-w$w.fastq"
  done
fi

if [ -n "${BBDUK:-}" ]; then
  for h in 1 2; do
    tools="$tools bbduk-hdist$h"
    common=(ref="$primers" k=12 mink=8 hdist="$h" rcomp=t threads=4 -Xmx2g
            ow=t)
    "$BBDUK" in="$reads" out="$out/bbduk-left$h.fastq.gz" ktrim=l \
      restrictleft=100 "${common[@]}" 2> "$out/bbduk-left$h.log"
    "$BBDUK" in="$out/bbduk-left$h.fastq.gz" \
      out="$out/bbduk-hdist$h.fastq.gz" ktrim=r restrictright=100 \
      minlen=1 "${common[@]}" 2> "$out/bbduk-right$h.log"
  done
fi

if [ -n "${PORECHOP:-}" ]; then
  adapters=$out/porechop-adapters.txt
  printf 'uni12A\n%s\n%s\nuni12G\n%s\n%s\n' "$uni12a" "$rc_mbt13" \
    "$uni12g" "$rc_mbt13" > "$adapters"
  tools="$tools porechop-abi porechop-abi-tuned"
  "$PORECHOP" -i "$reads" -o "$out/porechop-abi.fastq.gz" -cap "$adapters" \
    -ddb -t 4 > "$out/porechop-abi.log" 2>&1
  "$PORECHOP" -i "$reads" -o "$out/porechop-abi-tuned.fastq.gz" \
    -cap "$adapters" -ddb -t 4 --extra_end_trim 0 --end_threshold 60 \
    > "$out/porechop-abi-tuned.log" 2>&1
fi

if [ -n "${DORADO:-}" ]; then
  tools="$tools dorado"
  dp=$out/dorado-primers.fasta
  : > "$dp"
  for form in A:$uni12a G:$uni12g; do
    printf '>uni12%s_front\tet:Z:primer\tsk:Z:ANY\n%s\n' "${form%%:*}" \
      "${form#*:}" >> "$dp"
    printf '>uni12%s_rear\tet:Z:primer\tsk:Z:ANY\n%s\n' "${form%%:*}" \
      "$mbt13" >> "$dp"
  done
  input=$reads
  if [ $barcoded = 1 ]; then
    rm -rf "$out/dorado-demux"
    "$DORADO" demux --kit-name SQK-NBD114-24 --emit-fastq \
      -o "$out/dorado-demux" -t 4 "$reads" 2> "$out/dorado-demux.log"
    input=$out/dorado-demux.fastq
    find "$out/dorado-demux" -name '*.fastq' -exec cat {} + > "$input"
  fi
  "$DORADO" trim -k SQK-LSK114 --primer-sequences "$dp" --emit-fastq -t 4 \
    "$input" 2> "$out/dorado.log" | gzip > "$out/dorado.fastq.gz"
fi

if [ -n "${SAMTOOLS:-}" ] && [ -n "${MINIMAP2:-}" ]; then
  tools="$tools ampliconclip"
  genome=$out/genome.fasta
  cp "$sim/genome.fasta" "$genome"
  "$SAMTOOLS" faidx "$genome"
  awk -v OFS='\t' '{ print $1, 0, 12, "uni12_" $1, 0, "+"
                     print $1, $2 - 13, $2, "uni13_" $1, 0, "-" }' \
    "$genome.fai" > "$out/primers.bed"
  "$MINIMAP2" -ax map-ont -t 4 "$genome" "$reads" 2> "$out/minimap2.log" \
    | "$SAMTOOLS" view -b -F 0x904 - \
    | "$SAMTOOLS" sort -o "$out/mapped.bam" -
  "$SAMTOOLS" index "$out/mapped.bam"
  "$SAMTOOLS" ampliconclip --hard-clip --both-ends -b "$out/primers.bed" \
    -o "$out/clipped.bam" "$out/mapped.bam" 2> "$out/ampliconclip.log"
  # the soft-clipped ends (outside the alignment) are removed as well
  "$SAMTOOLS" view "$out/clipped.bam" | awk -F'\t' -v OFS='\n' '
    { c = $6; s = $10; q = $11; l = 0; r = 0
      if (match(c, /^[0-9]+S/)) l = substr(c, RSTART, RLENGTH - 1) + 0
      if (match(c, /[0-9]+S$/)) r = substr(c, RSTART, RLENGTH - 1) + 0
      n = length(s) - l - r
      if (n > 0)
        print "@" $1, substr(s, l + 1, n), "+", substr(q, l + 1, n) }' \
    | gzip > "$out/ampliconclip.fastq.gz"
fi

if [ -n "${BARBELL:-}" ] && [ $barcoded = 1 ]; then
  tools="$tools barbell"
  bb=$out/barbell
  rm -rf "$bb"
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
  "$BARBELL" annotate -q "$bb/left.fasta" "$bb/right.fasta" -b Ftag Rtag \
    -i "$reads" -o "$bb/anno.tsv" -t 4 > "$bb/annotate.log" 2>&1
  "$BARBELL" filter -i "$bb/anno.tsv" -f "$bb/filters.txt" \
    -o "$bb/filtered.tsv" > "$bb/filter.log" 2>&1
  "$BARBELL" trim -i "$bb/filtered.tsv" -r "$reads" -o "$bb/trimmed" \
    --no-orientation --no-flanks --sort-labels > "$bb/trim.log" 2>&1
  cat "$bb"/trimmed/*.fastq | gzip > "$out/barbell.fastq.gz"
fi

score="python3 -I $here/score_reads.py"
compare="python3 -I $here/compare_consensus.py"
header=
: > "$out/reads.tsv"
: > "$out/consensus.tsv"
printf 'run\tcodon\n' > "$out/site.tsv"
for tool in $tools; do
  $score "$sim" "$out/$tool.fastq.gz" "$tool" $header >> "$out/reads.tsv"
  for mode in plain elong; do
    run=$tool
    opts=()
    if [ "$mode" = elong ]; then
      run=${tool}_elong
      opts=(-c "$out/elong.sh")
    fi
    rm -rf "$out/irma/$run"
    mkdir -p "$out/irma"
    "$IRMA" FLU-minion "$out/$tool.fastq.gz" "$out/irma/$run" "${opts[@]}" \
      > "$out/irma/$run.log" 2>&1
    $compare "$sim" "$out/irma/$run" "$run" $header >> "$out/consensus.tsv"
    header=--no-header
    ha=$(grep -v '>' "$out/irma/$run"/A_HA_*.fasta | tr -d '\n')
    codon=$(echo "$ha" | { grep -o 'GGAGC[A-Z]\{4\}GCAGGG' || true; } \
            | head -1 | cut -c7-9)
    printf '%s\t%s\n' "$run" "${codon:-not found}" >> "$out/site.tsv"
  done
done
cat "$out/reads.tsv" "$out/site.tsv"
# the consensus, summed over the segments
awk -F'\t' 'NR > 1 { m5[$1] += $5; m3[$1] += $6; x5[$1] += $7;
  x3[$1] += $8; e[$1] += $9 + $10 + $11; a[$1] += $12; order[NR] = $1 }
  END { print "run\tmissing_5\tmissing_3\textra_5\textra_3\terrors\tambiguous"
        for (i = 2; i <= NR; i++) if (!(order[i] in done)) {
          r = order[i]; done[r] = 1
          print r "\t" m5[r] "\t" m3[r] "\t" x5[r] "\t" x3[r] "\t" e[r] \
            "\t" a[r] } }' "$out/consensus.tsv"
