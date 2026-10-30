#!/bin/bash
#
# Download the genome of the simulations from NCBI (E-utilities) and
# write it, with the date of the download, to <outdir>:
#
# - segments 1-3 and 5-8: A/Puerto Rico/8/1934 (H1N1), RefSeq
#   NC_002016.1-NC_002023.1 (without NC_002017.1, its HA), complete
#   sequences with both conserved termini;
# - segment 4: the HA of A/chicken/Cambodia/X0815302/2013 (H5N1),
#   KF918483.1, complete with both termini, which carries the Uni12-like
#   site of the I348K artefact (see benchmark/README.md).
#
# The combination is constructed for the simulation (PR8 internal genes
# with an H5 HA); it is not a sequenced isolate.
#
# Usage: benchmark/make-genome.sh <outdir>
#

set -euo pipefail
outdir=$1
mkdir -p "$outdir"
eutils=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
pr8=NC_002016.1,NC_002018.1,NC_002019.1,NC_002020.1,NC_002021.1
pr8=$pr8,NC_002022.1,NC_002023.1
ha=KF918483.1
query="db=nuccore&rettype=fasta&retmode=text"
curl -sSf "$eutils?$query&id=$pr8" > "$outdir/pr8.fasta"
curl -sSf "$eutils?$query&id=$ha" \
  | sed '1s/$/ [used as segment 4]/' > "$outdir/ha.fasta"
cat "$outdir/pr8.fasta" "$outdir/ha.fasta" | grep -v '^$' \
  > "$outdir/genome.fasta"
rm "$outdir/pr8.fasta" "$outdir/ha.fasta"
date -u +%FT%TZ > "$outdir/fetch_date.txt"
grep '>' "$outdir/genome.fasta"
