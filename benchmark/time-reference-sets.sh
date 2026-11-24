#!/bin/bash
# (c) 2026 Giorgio Gonnella
# Wall time and peak memory of flucut with each reference set, on a small
# and a large file of reads. A second flucut binary can be given to
# compare two builds on the same inputs (e.g. to check that a change to
# the voting rule costs no time).
#
# usage: time-reference-sets.sh <reads dir> <flucut> [other flucut]
#
# <reads dir> holds big.fastq.gz and small.fastq.gz. The reference sets
# are those of the cache ($XDG_CACHE_HOME/flucut or ~/.cache/flucut):
# the built-in sets and, if present, the set built from IRMA's seeds by
# scripts/irma_references.py. Writes a table of one row per run to
# standard output.
set -u

if [ $# -lt 2 ]; then
  sed -n '3,14p' "$0" | sed 's/^# \?//'
  exit 2
fi
reads_dir=$1
new=$2
old=${3:-}

cache=${XDG_CACHE_HOME:-$HOME/.cache}/flucut
sets=()
for spec in "v2:$cache/flu-a-references-v2.fasta" \
            "v3:$cache/flu-a-references-v3.fasta" \
            "irma:$cache/irma-flu-references.fasta"; do
  [ -f "${spec#*:}" ] && sets+=("$spec")
done
if [ ${#sets[@]} -eq 0 ]; then
  echo "no reference set in $cache" >&2
  exit 1
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

echo -e "binary\treads\treference_set\trepeat\tseconds\tKB"
for reads in big small; do
  [ -f "$reads_dir/$reads.fastq.gz" ] || continue
  for spec in "${sets[@]}"; do
    set_name=${spec%%:*}
    refs=${spec#*:}
    for bin_spec in "new:$new" "old:$old"; do
      label=${bin_spec%%:*}
      prog=${bin_spec#*:}
      [ -n "$prog" ] && [ -x "$prog" ] || continue
      for rep in 1 2 3; do
        t=$( { /usr/bin/time -f "%e %M" "$prog" -j 8 -r "$refs" \
               -o "$work/out" "$reads_dir/$reads.fastq.gz" \
               > /dev/null; } 2>&1 | tail -1 )
        echo -e "$label\t$reads\t$set_name\t$rep\t${t% *}\t${t#* }"
        rm -rf "$work/out"
      done
    done
  done
done
