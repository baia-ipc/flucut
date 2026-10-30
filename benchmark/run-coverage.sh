#!/bin/bash
#
# The tools at low coverage: reads simulated as in run-irma-checks.sh,
# with 10, 20 and 40 reads per segment (three seeds each), trimmed by
# the tools of run-tools.sh (same commands), assembled by IRMA
# FLU-minion (as it is); for each run, summed over the segments: the
# segments assembled, the differences of the majority consensus from the
# truth (substitutions, insertions, deletions; compare_consensus.py
# --plain), its ambiguous bases, its differences at the internal priming
# site, the N of the amended consensus (where fewer than 50 reads support
# a base, MIN_CONS_SUPPORT of FLU-minion), and the codon of the I348K
# site in the majority consensus of HA.
#
# Usage: benchmark/run-coverage.sh <genome.fasta> <workdir> <flucut> \
#          <IRMA> <irma-core>   (cutadapt in the PATH)
# Output: <workdir>/coverage.tsv, also shown.
#

set -euo pipefail
genome=$1
work=$2
flucut=$3
irma=$4
irma_core=$5
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$work"
primers=$work/irma-core-primers.fasta
printf '>uni12A\nAGCAAAAGCAGG\n>uni12G\nAGCGAAAGCAGG\n' > "$primers"
printf '>mbtuni13\nACGCGTGATCAGTAGAAACAAGG\n' >> "$primers"
table=$work/coverage.tsv
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' reads_per_segment seed tool \
  segments errors ambiguous site_errors amended_n codon > "$table"
for n in 10 20 40; do
  for seed in 1 2 3; do
    run=$work/n${n}_s$seed
    mkdir -p "$run"
    python3 -I "$here/simulate.py" "$genome" "$run/sim" \
      --reads-per-segment "$n" --seed "$seed"
    reads=$run/sim/reads.fastq.gz
    cp "$reads" "$run/untrimmed.fastq.gz"
    "$flucut" -q -o "$run/flucut" "$reads"
    mv "$run/flucut/reads.fastq.gz" "$run/flucut.fastq.gz"
    cutadapt -j 1 -e 0.2 --revcomp -n 2 -m 1 -g AGCRAAAGCAGG \
      -a CCTTGTTTCTACTGATCACGCGT -o "$run/cutadapt.fastq.gz" "$reads" \
      > "$run/cutadapt.log"
    "$irma_core" trimmer "$reads" -o "$run/irma-core-w70.fastq" \
      -P "$primers" --p-kmer-length 12 --p-fuzzy --p-restrict 70 -n 1 \
      2> "$run/irma-core-w70.log"
    gzip -f "$run/irma-core-w70.fastq"
    for tool in untrimmed flucut cutadapt irma-core-w70; do
      rm -rf "$run/irma/$tool"
      mkdir -p "$run/irma"
      "$irma" FLU-minion "$run/$tool.fastq.gz" "$run/irma/$tool" \
        > "$run/irma/$tool.log" 2>&1 || true
      stats=$(python3 -I "$here/compare_consensus.py" --plain --no-header \
                "$run/sim" "$run/irma/$tool" "$tool" | awk -F'\t' '
        $3 != "NA" { s++; e += $9 + $10 + $11; a += $12 }
        $13 != "NA" { site += $13 }
        END { printf "%d\t%d\t%d\t%d", s, e, a, site }')
      amended_n=$(python3 -I "$here/compare_consensus.py" --no-header \
                    "$run/sim" "$run/irma/$tool" "$tool" | awk -F'\t' '
        $3 != "NA" { n += $14 } END { printf "%d", n }')
      codon=not_assembled
      if ls "$run/irma/$tool"/A_HA_*.fasta > /dev/null 2>&1; then
        ha=$(grep -v '>' "$run/irma/$tool"/A_HA_*.fasta | tr -d '\n')
        codon=$(echo "$ha" | { grep -o 'GGAGC[A-Z]\{4\}GCAGGG' || true; } \
                | head -1 | cut -c7-9)
        codon=${codon:-not_found}
      fi
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$n" "$seed" "$tool" "$stats" \
        "$amended_n" "$codon" >> "$table"
    done
  done
done
cat "$table"
