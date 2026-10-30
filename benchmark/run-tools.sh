#!/bin/bash
#
# The comparison of the tools on the simulated reads of
# run-irma-checks.sh (<workdir>/sim): each tool trims the reads, the
# trimmed reads are scored against the truth (score_reads.py), and IRMA
# FLU-minion assembles them, as it is and with the reference elongation
# (SKIP_E=0); the consensus is compared with the truth
# (compare_consensus.py), and the codon of the I348K site is read in the
# plain consensus of HA.
#
# The commands of the tools are those of benchmark/RESULTS.md: cutadapt
# with the two primers non-anchored (-g, -a), --revcomp, -e 0.2; the
# IRMA-core trimmer with the primers as k-mers (k 12, the length of
# Uni12; one mismatch), windows of 40 (the default of MIRA-NF for custom
# primers) and 70 bases.
#
# Usage: benchmark/run-tools.sh <workdir> <flucut> <IRMA> <irma-core>
#   (cutadapt in the PATH)
# Output: <workdir>/tools/{reads,consensus,site}.tsv, also shown.
#

set -euo pipefail
work=$1
flucut=$2
irma=$3
irma_core=$4
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
sim=$work/sim
reads=$sim/reads.fastq.gz
out=$work/tools
mkdir -p "$out"
printf 'SKIP_E=0\n' > "$out/elong.sh"
primers=$out/irma-core-primers.fasta
printf '>uni12A\nAGCAAAAGCAGG\n>uni12G\nAGCGAAAGCAGG\n' > "$primers"
printf '>mbtuni13\nACGCGTGATCAGTAGAAACAAGG\n' >> "$primers"

# the trimmed reads of each tool, in $out/<tool>.fastq.gz
cp "$reads" "$out/untrimmed.fastq.gz"
"$flucut" -q -o "$out/flucut" "$reads"
mv "$out/flucut/reads.fastq.gz" "$out/flucut.fastq.gz"
cutadapt -j 1 -e 0.2 --revcomp -n 2 -m 1 -g AGCRAAAGCAGG \
  -a CCTTGTTTCTACTGATCACGCGT -o "$out/cutadapt.fastq.gz" "$reads" \
  > "$out/cutadapt.log"
for w in 40 70; do
  "$irma_core" trimmer "$reads" -o "$out/irma-core-w$w.fastq" \
    -P "$primers" --p-kmer-length 12 --p-fuzzy --p-restrict "$w" -n 1 \
    -v 2> "$out/irma-core-w$w.log"
  gzip -f "$out/irma-core-w$w.fastq"
done

tools="untrimmed flucut cutadapt irma-core-w40 irma-core-w70"
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
    "$irma" FLU-minion "$out/$tool.fastq.gz" "$out/irma/$run" "${opts[@]}" \
      > "$out/irma/$run.log" 2>&1
    $compare "$sim" "$out/irma/$run" "$run" $header >> "$out/consensus.tsv"
    header=--no-header
    ha=$(grep -v '>' "$out/irma/$run"/A_HA_*.fasta | tr -d '\n')
    codon=$(echo "$ha" | grep -o 'GGAGCT[ACGT]\{3\}GCAGGG' | head -1 \
            | cut -c7-9)
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
