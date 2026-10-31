#!/bin/bash
#
# Run a program test case of flucut and write a digest of its results to
# the standard output, compared by the Makefile (test_programs) with
# testdata/expected/<case>.txt.
#
# A case is a line of testdata/cases.tsv: its name, a tab and the
# arguments of flucut (the output directory is added). The digest has
# the exit status, the messages without their time stamps, and, for each
# output file, the md5 of its decompressed content and, for each read,
# its title, its length and its first and last 15 bases; then the
# summary table, the lengths table and the JSON summary (with the time
# of writing replaced by DATE); then, for the HTML report, its last line
# and its headings (its content is checked through the JSON summary,
# which has the same content).
#
# Usage: tests/run-case.sh <flucut> <case> "<valgrind command or empty>"
#

set -euo pipefail
export LC_ALL=C
# the built-in reference set is a local file in the tests (no download)
export FLUCUT_BUILTIN_REFERENCES=testdata/ref-ha.fasta
prog=$1
name=$2
valgrind=${3:-}
args=$(awk -F'\t' -v c="$name" '$1 == c { print $2 }' testdata/cases.tsv)
outdir=$(mktemp -d "${TMPDIR:-/tmp}/flucut-case-XXXXXX")
trap 'rm -rf "$outdir"' EXIT
status=0
# (the arguments are split at the spaces: no path of a case has one)
# shellcheck disable=SC2086
$valgrind "$prog" -o "$outdir" $args 2> "$outdir/messages" || status=$?
echo "exit status: $status"
echo "== messages =="
sed -e 's/^# [0-9-]* [0-9:]* | //' -e "s|$outdir|OUTDIR|g" \
  "$outdir/messages"
for f in "$outdir"/*.fastq.gz; do
  [ -e "$f" ] || continue
  echo "== $(basename "$f") =="
  echo "md5: $(gzip -dc "$f" | md5sum | cut -d' ' -f1)"
  gzip -dc "$f" | awk 'NR % 4 == 1 { title = $0 }
    NR % 4 == 2 { n = length($0); print title
                  print n, substr($0, 1, 15) "..." substr($0, n - 14) }'
done
if [ -e "$outdir/summary.tsv" ]; then
  echo "== summary.tsv =="
  cat "$outdir/summary.tsv"
fi
if [ -e "$outdir/lengths.tsv" ]; then
  echo "== lengths.tsv =="
  cat "$outdir/lengths.tsv"
fi
if [ -e "$outdir/flucut_summary.json" ]; then
  echo "== flucut_summary.json =="
  sed -e 's/^\(  "updated": \).*/\1"DATE",/' -e "s|$outdir|OUTDIR|g" \
    "$outdir/flucut_summary.json"
fi
html="$outdir/flucut_report.html"
if [ -e "$html" ]; then
  echo "== flucut_report.html =="
  tail -n 1 "$html"
  grep -oE '<h[23][^>]*>[^<]*</h[23]>' "$html"
fi
