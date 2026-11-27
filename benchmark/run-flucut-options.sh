#!/bin/bash
# (c) 2026 Giorgio Gonnella
#
# The read-level scores of flucut alone (score_reads.py) on simulated
# reads (a directory written by simulate.py), by its options:
# - the references: default (the built-in set), -r off, -r off --virus a
#   (the primers of influenza A alone), -r with the simulated genome;
# - the adapter search: each of them also with --adapters none;
# - the errors allowed in a primer: -e 0.1 to 0.35 instead of auto, with
#   -r off and with -r off --virus a.
# The error rate measured by -e auto and the value it chose are in the
# log of each run (<output dir>/<label>.log).
#
# Usage: benchmark/run-flucut-options.sh <simulation dir> <output dir>
# FLUCUT is given in the environment.
# Output: <output dir>/flucut_options.tsv, also shown.

set -euo pipefail
sim=$1
out=$2
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
reads=$sim/reads.fastq.gz
mkdir -p "$out"
table=$out/flucut_options.tsv
header=
: > "$table"
run() {  # <label> <flucut options>...
  local label=$1
  shift
  "$FLUCUT" -v -L "$out/$label.log" -o "$out/$label" "$@" "$reads"
  python3 -I "$here/score_reads.py" "$sim" "$out/$label/reads.fastq.gz" \
    "$label" $header >> "$table"
  header=--no-header
}
run default
run r-off -r off
run r-off-a -r off --virus a
run genome -r "$sim/genome.fasta"
run default-no-adapters --adapters none
run r-off-no-adapters -r off --adapters none
run r-off-a-no-adapters -r off --virus a --adapters none
run genome-no-adapters -r "$sim/genome.fasta" --adapters none
for e in 0.1 0.15 0.2 0.25 0.3 0.35; do
  run r-off-e$e -r off -e $e
  run r-off-a-e$e -r off --virus a -e $e
done
cat "$table"
