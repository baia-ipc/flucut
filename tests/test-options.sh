#!/bin/bash
#
# Validation of the command line of flucut: each invalid command line
# must be rejected with exit status 1 (the Makefile gives valgrind
# another exit status for its own errors); the help and the version must
# succeed.
#
# Usage: tests/test-options.sh <flucut> "<valgrind command or empty>"
#

set -uo pipefail
prog=$1
valgrind=${2:-}
outdir=$(mktemp -d "${TMPDIR:-/tmp}/flucut-options-XXXXXX")
trap 'rm -rf "$outdir"' EXIT
reads=testdata/reads.fastq
refs=testdata/ref-ha.fasta
failed=0

# run flucut with the arguments after the first; the first is the exit
# status expected
expect() {
  local expected=$1
  shift
  local status=0
  # shellcheck disable=SC2086
  $valgrind "$prog" -o "$outdir/out" "$@" > /dev/null 2>&1 || status=$?
  if [ "$status" -ne "$expected" ]; then
    echo "[FAIL] options test: flucut $* (exit status $status," \
      "expected $expected)"
    failed=1
  fi
}

expect 0 -h
expect 0 -V
# no input
expect 1
# keeping the termini needs the references
expect 1 -T terminus "$reads"
expect 1 -T both -r "$refs" "$reads"
expect 1 -s no-such-scheme "$reads"
expect 1 -e -1 "$reads"
expect 1 -e 1.5 "$reads"
expect 1 -e x "$reads"
expect 1 -e autox "$reads"
expect 1 -w -1 "$reads"
expect 1 -m 0 "$reads"
expect 1 -j 0 "$reads"
expect 1 --adapters nb,xx "$reads"
expect 1 --adapters "" "$reads"
expect 1 -j 1000 "$reads"
expect 1 --probe-len 10 -r "$refs" "$reads"
expect 1 --probe-errors 1 -r "$refs" "$reads"
expect 1 --terminal-window -5 -r "$refs" "$reads"
expect 1 --no-such-option "$reads"
# inputs: missing, not a file of reads, two with the same output
expect 1 testdata/no-such-file.fastq
expect 1 testdata/cases.tsv
expect 1 "$reads" "$reads"
# two groups with the same name (barcode01)
expect 1 testdata/barcode01 testdata/fastq_pass
# references: missing, not FASTA
expect 1 -r testdata/no-such-file.fasta "$reads"
expect 1 -r "$reads" "$reads"
expect 1 -o "" "$reads"

if [ "$failed" -eq 0 ]; then
  echo "[OK] options tests"
fi
exit "$failed"
