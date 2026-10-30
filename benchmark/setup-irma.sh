#!/bin/bash
#
# Install IRMA from its sources into <dir>, as done for the benchmark:
# the IRMA repository, LABEL (IRMA uses its align2model binary and its
# scripts), and IRMA-core and IRMA-viz built with Rust nightly (their
# release archives are the usual way; they were not reachable from the
# environment of the benchmark). The versions are those pinned by IRMA
# (.package-core.sh, .package-label.sh, .package-viz.sh of IRMA).
#
# Requires: git, cargo with rustup (a nightly toolchain is installed).
#
# Usage: benchmark/setup-irma.sh <dir>
# then:  <dir>/irma/IRMA FLU-minion reads.fastq.gz <sample>
#

set -euo pipefail
dir=$1
irma_commit=1311c81
core=v0.11.0
viz=v0.2.2
label=v0.7.2
mkdir -p "$dir"
cd "$dir"
gh=https://github.com/CDCgov
git clone -q "$gh/irma" irma
git -C irma checkout -q "$irma_commit"
git clone -q --depth 1 --branch "$label" "$gh/label" label
git clone -q --depth 1 --branch "$core" "$gh/irma-core" irma-core
git clone -q --depth 1 --branch "$viz" "$gh/irma-viz" irma-viz
rustup toolchain install nightly --profile minimal
(cd irma-core && cargo +nightly build --profile prod)
(cd irma-viz && cargo +nightly build --profile prod)
suffix="$(uname -s)_$(uname -m)"
cp irma-core/target/prod/irma-core "irma/IRMA_RES/scripts/irma-core_$suffix"
cp irma-viz/target/prod/irma-viz "irma/IRMA_RES/scripts/irma-viz_$suffix"
irma/IRMA_RES/scripts/irma-core_"$suffix" --version
