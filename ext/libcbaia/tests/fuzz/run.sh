#!/bin/bash
#
# Run a fuzzing harness (built by tests/fuzz/build.sh).
#
# Usage: tests/fuzz/run.sh <harness> <seconds> [<jobs>]
#   harness: reads, bam or whitelist
#
# The corpus is kept in obj/fuzz/corpus-<harness> (initialized with the
# seeds in tests/fuzz/seeds/<harness>); crashing inputs are written to
# obj/fuzz/crashes-<harness>/.
#

set -euo pipefail
HARNESS=$1
SECONDS_MAX=$2
JOBS=${3:-1}
CORPUS=obj/fuzz/corpus-$HARNESS
CRASHES=obj/fuzz/crashes-$HARNESS
mkdir -p $CORPUS $CRASHES
cp -n tests/fuzz/seeds/$HARNESS/* $CORPUS/ 2>/dev/null || true
# leak detection is disabled, as runs interrupted by exit() leave memory
# allocated (see fuzz-common.h); for the same reason, the fork mode is used
# (the inputs are run in short-lived child processes); the output of the
# target is suppressed
ASAN_OPTIONS=detect_leaks=0 bin/fuzz/fuzz-$HARNESS $CORPUS \
  -max_total_time=$SECONDS_MAX -detect_leaks=0 -close_fd_mask=3 \
  -max_len=8192 -timeout=10 -artifact_prefix=$CRASHES/ -fork=$JOBS \
  -print_final_stats=1 -rss_limit_mb=4096
