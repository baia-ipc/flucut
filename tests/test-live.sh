#!/bin/bash
# (c) 2026 Giorgio Gonnella
#
# Test of the live mode (--live; from the test of the real-time mode of
# Lumipore, tests/test-watch.sh): files are written into the input
# directory while flucut is running, as by MinKNOW (including a file
# written in two parts, and a barcode directory created later); the
# outputs and the tables must be the same as those of a normal run on the
# final input. Then: --resume with a new file, restart without --resume,
# interruption with SIGINT, end of the run (final summary of MinKNOW),
# -e auto with no file at the start, a file created empty.
#
# Usage: tests/test-live.sh <flucut command> [<prefix command>]
#

set -euo pipefail
export LC_ALL=C
FLUCUT=$1
PREFIX=${2:-}
SRC=testdata/fastq_pass
TMP=$(mktemp -d)
# at the exit (also after a failure): the background flucut stopped
trap 'kill $(jobs -p) 2> /dev/null || true; rm -rf $TMP' EXIT
RUN=$TMP/fastq_pass

fail() { echo "FAIL: $*" >&2; exit 1; }

# wait until a condition holds: <command> is run every 0.1 s, for at most
# WAIT_TIMEOUT seconds (instead of fixed pauses, which can be too short,
# e.g. under valgrind, or longer than needed); usage:
# wait_for <description> <command> [<arguments>]
WAIT_TIMEOUT=120
wait_for() {
  local what=$1 i
  shift
  for (( i = 0; i < WAIT_TIMEOUT * 10; i++ )); do
    "$@" && return 0
    sleep 0.1
  done
  fail "timeout waiting for: $what"
}

# conditions for wait_for: a file contains a string; a file contains at
# least <n> lines matching a pattern
contains() { grep -q "$2" "$1" 2> /dev/null; }
n_lines() { [ "$(grep -c "$2" "$1" 2> /dev/null)" -ge "$3" ]; }

# the files processed, in the messages
n_processed() { grep -c '| barcode0[12]: [^ ]*\.fastq' "$1" || true; }

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

# write a file in two parts, with a pause (as a file being written); the
# pause is longer than --min-age (1 s) of the test, so that the incomplete
# file is old enough to be processed and only its content holds it back
write_slowly() {
  local src=$1 dst=$2 half=$(( $(wc -c < $1) / 2 ))
  head -c $half $src > $dst
  sleep 3
  tail -c +$(( half + 1 )) $src >> $dst
}

# 1. files written during the run, after the first scan (when the reports
#    of the run in progress are written); the scans are 6 s apart (the
#    --timeout); the pause of 2 s places the scan at 6 s while the second
#    file of barcode01 is written (old enough, but incomplete)
mkdir -p $RUN/barcode01
$PREFIX $FLUCUT -j 2 -o $TMP/out --live --timeout 6s --min-age 1 $RUN \
  2> $TMP/live.err &
pid=$!
wait_for "live: first scan" test -f $TMP/out/flucut_summary.json
grep -q '"running": true' $TMP/out/flucut_summary.json || \
  fail "summary of the run in progress not running"
cp $SRC/barcode02/part2.fastq.gz $RUN/barcode01/a.fastq.gz
sleep 2
write_slowly $SRC/barcode01/part1.fastq $RUN/barcode01/b.fastq
mkdir -p $RUN/barcode02
write_slowly testdata/reads.fastq $RUN/barcode02/a.fastq
wait $pid || fail "live mode exited with an error: $(cat $TMP/live.err)"
[ "$(n_processed $TMP/live.err)" -eq 3 ] || \
  fail "not 3 files processed: $(cat $TMP/live.err)"
grep -q 'modified after being processed' $TMP/live.err && \
  fail "a file was processed before being complete"
$FLUCUT -j 2 -o $TMP/dir $RUN 2> /dev/null
compare $TMP/out $TMP/dir
grep -q '"running": false' $TMP/out/flucut_summary.json || \
  fail "summary not final"
grep -q 'http-equiv="refresh"' $TMP/out/flucut_report.html && \
  fail "final report reloads itself"
echo "[OK] live: files written during the run"

# 2. resume, with a new file
cp testdata/reads.fastq $RUN/barcode01/c.fastq
$PREFIX $FLUCUT -o $TMP/out --live --timeout 1s --min-age 0 --resume $RUN \
  2> $TMP/resume.err
[ "$(n_processed $TMP/resume.err)" -eq 1 ] || \
  fail "resume: not only the new file processed: $(cat $TMP/resume.err)"
rm -rf $TMP/dir
$FLUCUT -o $TMP/dir $RUN 2> /dev/null
compare $TMP/out $TMP/dir
for g in barcode01 barcode02; do
  [ -s $TMP/out/$g.flucut-files ] && [ -s $TMP/out/$g.flucut-counts ] || \
    fail "resume: state files of $g missing"
done
[ -z "$(find $TMP/out -name '*.tmp')" ] || \
  fail "resume: temporary files left"
echo "[OK] live: resume"

# 3. without --resume, the outputs are written again from the beginning
$PREFIX $FLUCUT -o $TMP/out --live --timeout 1s --min-age 0 $RUN \
  2> /dev/null
compare $TMP/out $TMP/dir
echo "[OK] live: restart"

# 4. interruption (no time limit)
rm -rf $TMP/out
$PREFIX $FLUCUT -o $TMP/out --live --min-age 0 $RUN 2> $TMP/int.err &
pid=$!
wait_for "interruption: the 4 files processed" \
  n_lines $TMP/int.err '| barcode0[12]: ' 4
kill -INT $pid
wait $pid || fail "interrupted live mode exited with an error"
grep -q 'Interrupted' $TMP/int.err || \
  fail "interrupted: not reported: $(cat $TMP/int.err)"
compare $TMP/out $TMP/dir
echo "[OK] live: interruption"

# 5. end of the run: when MinKNOW writes the final summary in the run
#    directory (the parent of fastq_pass), the files are processed without
#    waiting for their minimum age, then flucut stops (no time limit); the
#    run ends after the first scan
rm -rf $TMP/end $TMP/end_out
mkdir -p $TMP/end/fastq_pass/barcode01
cp $SRC/barcode01/part1.fastq $TMP/end/fastq_pass/barcode01/a.fastq
$PREFIX $FLUCUT -j 2 -o $TMP/end_out --live --min-age 3600 \
  $TMP/end/fastq_pass 2> $TMP/end.err &
pid=$!
wait_for "end of the run: first scan" test -f $TMP/end_out/flucut_summary.json
cp $SRC/barcode02/part2.fastq.gz $TMP/end/fastq_pass/barcode01/b.fastq.gz
touch $TMP/end/final_summary_FAX12345_abcdef12_0123abcd.txt
wait $pid || fail "end of the run: error: $(cat $TMP/end.err)"
[ "$(n_processed $TMP/end.err)" -eq 2 ] || \
  fail "end of the run: not 2 files processed: $(cat $TMP/end.err)"
grep -q 'End of the run' $TMP/end.err || \
  fail "end of the run: not reported: $(cat $TMP/end.err)"
echo "[OK] live: end of the run"

# 6. -e auto, with no file at the start: decided at the first file, as
#    in a normal run
rm -rf $TMP/auto $TMP/auto_out $TMP/auto_dir
mkdir -p $TMP/auto
$PREFIX $FLUCUT -v -o $TMP/auto_out --live --timeout 5s --min-age 0 \
  $TMP/auto 2> $TMP/auto.err &
pid=$!
wait_for "auto: start without files" \
  test -f $TMP/auto_out/flucut_summary.json
grep -q '"max_errors": {"mode": "auto", "value": null' \
  $TMP/auto_out/flucut_summary.json || fail "auto: decided without a file"
cp testdata/reads.fastq $TMP/auto/
wait $pid || fail "auto: error: $(cat $TMP/auto.err)"
grep -q 'errors in a primer: 0.2 of its length (-e auto)' $TMP/auto.err || \
  fail "auto: not decided at the first file: $(cat $TMP/auto.err)"
grep -q '"max_errors": {"mode": "auto", "value": 0.2, ' \
  $TMP/auto_out/flucut_summary.json || \
  fail "auto: decision not in the summary"
$FLUCUT -o $TMP/auto_dir $TMP/auto 2> /dev/null
compare $TMP/auto_out $TMP/auto_dir
echo "[OK] live: -e auto"

# 7. a file created empty and written later (as by a copy): not processed
#    while empty, then processed completely
rm -rf $TMP/empty $TMP/empty_out
mkdir -p $TMP/empty
$PREFIX $FLUCUT -o $TMP/empty_out --live --timeout 4s --min-age 0 \
  $TMP/empty 2> $TMP/empty.err &
pid=$!
: > $TMP/empty/reads.fastq
wait_for "empty file: first scan" test -f $TMP/empty_out/flucut_summary.json
sleep 0.5
cat testdata/reads.fastq > $TMP/empty/reads.fastq
wait $pid || fail "empty file: error: $(cat $TMP/empty.err)"
grep -q '"in": 12,' $TMP/empty_out/flucut_summary.json || \
  fail "empty file: not all reads processed: $(cat $TMP/empty.err)"
echo "[OK] live: file created empty"
