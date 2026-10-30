#!/bin/bash
#
# The IRMA checks of the benchmark (see benchmark/RESULTS.md): on the
# simulated reads, the IRMA consensus of
#   raw        the untrimmed reads
#   flucut     flucut with the default trimming (every primer removed)
#   notags     the same, with --no-tags
#   terminus   flucut --trim terminus with the genome as references
# each with IRMA FLU-minion as it is, and with the reference elongation
# on (SKIP_E=0, suffix _elong), compared with the truth
# (benchmark/compare_consensus.py). Also the codon of the I348K site in
# the plain (majority) consensus of HA.
#
# Usage: benchmark/run-irma-checks.sh <genome.fasta> <workdir> <flucut> \
#          <IRMA>
# Output: <workdir>/consensus.tsv and <workdir>/site.tsv, also shown.
#

set -euo pipefail
genome=$1
work=$2
flucut=$3
irma=$4
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$work"
python3 -I "$here/simulate.py" "$genome" "$work/sim"
reads=$work/sim/reads.fastq.gz
"$flucut" -q -o "$work/flucut" "$reads"
"$flucut" -q --no-tags -o "$work/notags" "$reads"
"$flucut" -q -T terminus -r "$work/sim/genome.fasta" -o "$work/terminus" \
  "$reads"
printf 'SKIP_E=0\n' > "$work/elong.sh"
python_compare="python3 -I $here/compare_consensus.py"
first=1
: > "$work/consensus.tsv"
printf 'run\thas_codon\n' > "$work/site.tsv"
for name in raw flucut notags terminus; do
  if [ "$name" = raw ]; then input=$reads; else
    input=$work/$name/reads.fastq.gz; fi
  for mode in plain elong; do
    run=$name; opts=()
    if [ "$mode" = elong ]; then run=${name}_elong; opts=(-c "$work/elong.sh")
    fi
    rm -rf "$work/irma/$run"
    mkdir -p "$work/irma"
    "$irma" FLU-minion "$input" "$work/irma/$run" "${opts[@]}" \
      > "$work/irma/$run.log" 2>&1
    if [ $first -eq 1 ]; then
      $python_compare "$work/sim" "$work/irma/$run" "$run" \
        >> "$work/consensus.tsv"
      first=0
    else
      $python_compare --no-header "$work/sim" "$work/irma/$run" "$run" \
        >> "$work/consensus.tsv"
    fi
    # the codon of the site (template ATA: Ile; primer AAA: Lys), found
    # by its context in the plain consensus of HA (the base before the
    # codon can be the R of the primer, A or G, instead of the T)
    ha=$(grep -v '>' "$work/irma/$run"/A_HA_*.fasta | tr -d '\n')
    codon=$(echo "$ha" | { grep -o 'GGAGC[A-Z]\{4\}GCAGGG' || true; } \
            | head -1 | cut -c7-9)
    printf '%s\t%s\n' "$run" "${codon:-not found}" >> "$work/site.tsv"
  done
done
cat "$work/consensus.tsv" "$work/site.tsv"
