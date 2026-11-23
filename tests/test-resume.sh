#!/bin/bash
# (c) 2026 Giorgio Gonnella
#
# Test of --resume in a normal run (not the live mode): a run of two
# barcodes is killed after the first one (its input is a named pipe, so
# the run waits there, with the first barcode complete); run again with
# --resume, it must skip the first barcode and give the same outputs and
# tables as a run without interruption. Then: --resume in a complete
# output directory trims nothing and only writes the reports again; a
# run without --resume trims everything again.
#
# Usage: tests/test-resume.sh <flucut command> [<prefix command>]
#

set -euo pipefail
export LC_ALL=C
# the built-in reference set is a local file in the tests (no download)
export FLUCUT_BUILTIN_REFERENCES=testdata/ref-ha.fasta
FLUCUT=$1
PREFIX=${2:-}
SRC=testdata/fastq_pass
TMP=$(mktemp -d)
# at the exit (also after a failure): the background flucut stopped
trap 'kill -9 $(jobs -p) 2> /dev/null || true; rm -rf $TMP' EXIT

fail() { echo "FAIL: $*" >&2; exit 1; }

# wait (at most 120 s, e.g. under valgrind) until a file exists
wait_file() {
  local i
  for (( i = 0; i < 1200; i++ )); do
    [ -e "$1" ] && return 0
    sleep 0.1
  done
  fail "timeout waiting for $1"
}

# compare the outputs (reads and tables) of two output directories
compare() {
  local a=$1 b=$2
  [ "$(cd $a && ls *.fastq.gz)" = "$(cd $b && ls *.fastq.gz)" ] || \
    fail "different output files in $a and $b"
  for f in $(cd $a && ls *.fastq.gz); do
    cmp -s <(gzip -dc $a/$f) <(gzip -dc $b/$f) || fail "$f differs"
  done
  for f in summary.tsv lengths.tsv; do
    cmp -s $a/$f $b/$f || fail "$f differs"
  done
}

# the state files of the groups (kept after a complete run)
STATE=.flucut-state
has_state() {
  for g in barcode01 barcode02; do
    [ -s $1/$STATE/$g.flucut-files ] && [ -s $1/$STATE/$g.flucut-counts ] \
      || fail "state files of $g missing in $1"
  done
  ! ls $1/$STATE/*.tmp > /dev/null 2>&1 || fail "temporary files left"
}

# the input: two barcodes, with the same paths in all the runs
RUN=$TMP/fastq_pass
mkdir -p $RUN/barcode01 $RUN/barcode02
cp $SRC/barcode01/part1.fastq $RUN/barcode01/
cp $SRC/barcode02/part2.fastq.gz $RUN/barcode02/

# 1. the run without interruption
$PREFIX $FLUCUT -o $TMP/whole $RUN 2> $TMP/whole.err || \
  fail "run: error: $(cat $TMP/whole.err)"
has_state $TMP/whole

# 2. a run killed after barcode01: barcode02 is a named pipe, which
#    blocks the run until it is killed
mv $RUN/barcode02/part2.fastq.gz $TMP/part2.fastq.gz
mkfifo $RUN/barcode02/part2.fastq.gz
$PREFIX $FLUCUT -o $TMP/out $RUN 2> $TMP/killed.err &
wait_file $TMP/out/$STATE/barcode01.flucut-files
sleep 0.5
kill -9 $! 2> /dev/null || true
wait $! 2> /dev/null || true
[ -e $TMP/out/$STATE/barcode01.flucut-counts ] || \
  fail "killed run: counts of barcode01 not saved"
[ ! -e $TMP/out/$STATE/barcode02.flucut-files ] || \
  fail "killed run: barcode02 saved as complete"
rm $RUN/barcode02/part2.fastq.gz
mv $TMP/part2.fastq.gz $RUN/barcode02/

# 3. resumed: barcode01 skipped, the same outputs as the whole run
$PREFIX $FLUCUT -o $TMP/out --resume $RUN 2> $TMP/resume.err || \
  fail "resume: error: $(cat $TMP/resume.err)"
grep -q "barcode01: trimmed by an earlier run, skipped" $TMP/resume.err || \
  fail "resume: barcode01 not skipped: $(cat $TMP/resume.err)"
! grep -q "barcode02: trimmed by an earlier run" $TMP/resume.err || \
  fail "resume: barcode02 skipped"
compare $TMP/whole $TMP/out
has_state $TMP/out
echo "[OK] resume: interrupted run continued"

# 4. --resume after a complete run: both barcodes skipped, the reads not
#    written again, the reports written again with the same content
cp $TMP/out/barcode01.fastq.gz $TMP/barcode01.before
rm $TMP/out/summary.tsv $TMP/out/flucut_report.html
$PREFIX $FLUCUT -o $TMP/out --resume $RUN 2> $TMP/again.err || \
  fail "resume of a complete run: error: $(cat $TMP/again.err)"
[ "$(grep -c "trimmed by an earlier run, skipped" $TMP/again.err)" -eq 2 ] \
  || fail "resume of a complete run: not all skipped: $(cat $TMP/again.err)"
cmp -s $TMP/barcode01.before $TMP/out/barcode01.fastq.gz || \
  fail "resume of a complete run: the reads were written again"
[ -s $TMP/out/flucut_report.html ] || \
  fail "resume of a complete run: report not written"
compare $TMP/whole $TMP/out
echo "[OK] resume: complete run, only the reports written again"

# 5. without --resume, everything is trimmed again
$PREFIX $FLUCUT -o $TMP/out $RUN 2> $TMP/fresh.err || \
  fail "run again: error: $(cat $TMP/fresh.err)"
! grep -q "trimmed by an earlier run" $TMP/fresh.err || \
  fail "run again without --resume: a barcode skipped"
compare $TMP/whole $TMP/out
echo "[OK] resume: without --resume, all trimmed again"
