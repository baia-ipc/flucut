# flucut-posthoc

Measures, after the fact and from the mapping of every read, what the
report of flucut predicts for a sample.

flucut is fast because it does not map the reads: it classifies each read
end from a probe of 32 bases after the primer, and builds its coverage
plot from that. `flucut-posthoc` takes the time flucut does not: it maps
every read base by base to complete segment references, and from the
mappings it measures the same things — what the reads are, the origin of
each amplicon end, the internal priming sites, the amplicon lengths, the
coverage by read origin — and writes them as an HTML report of the same
kind. With the outputs of flucut it compares the two, read by read; with
the results of IRMA it shows what IRMA did with each kind of read.

It is separate from flucut (Python, in this directory) and does not
change it.

## Install

    mamba env create -f environment.yml
    conda activate flucut-posthoc

This provides minimap2, samtools and the Python packages (edlib, numpy,
docopt, loguru, pytest).

## Quick start

The reads of one sample, as they were before flucut:

    ./flucut-posthoc -o posthoc/sample1 sample1.fastq.gz

The same, compared with flucut and with IRMA runs before and after
flucut:

    ./flucut-posthoc -o posthoc/sample1 \
        -i irma/sample1_IRMA,irma/sample1_flucut_IRMA \
        -f trimmed/sample1.fastq.gz \
        -s trimmed/flucut_summary.json -R trimmed/flucut_report.html \
        sample1.fastq.gz

Open `posthoc/sample1/posthoc_report.html`. For several samples, one
table with links to their reports:

    ./flucut-posthoc-overview -o posthoc posthoc/*/posthoc_summary.json

A sample of a million reads takes about two minutes with 32 threads
(`-j`).

## Inputs

| input | option | used for |
|---|---|---|
| reads before flucut (FASTQ) | argument | everything measured: their primers are the evidence |
| IRMA output directories | `-i` | read counts, consensus, depth, which reads IRMA assembled |
| reads written by flucut | `-f` | flucut's call for each read (its tags) |
| `flucut_summary.json` | `-s` | flucut's counts and amplicon lengths |
| `flucut_report.html` | `-R` | flucut's coverage plots |
| segment references | `-r` | default: the set flucut downloads (`~/.cache/flucut`) |

Only the reads are required. If you only have the reads written by
flucut, give those: the coverage, the positions and the deletions are
measured all the same, but the primers are gone, so the origin of the
ends and the trimming cannot be checked.

## Outputs

In the output directory:

- `posthoc_report.html`: the report (self-contained);
- `posthoc_summary.json`: the same content for programs;
- `reads.tsv.gz`: one line per amplicon of each read;
- `offtarget_amplicons.fasta`: a typical read of each group of
  amplicons which are not influenza.

## Documentation

- [Tutorial](docs/tutorial.md): one sample from start to finish, then a
  whole run.
- [Reading the report](docs/report.md): each part of the report, and the
  output files.
- [How it works](docs/method.md): the method, how it was checked, and
  its limits.

## Tests

    python -m pytest tests

The tests include the whole analysis on simulated reads with a known
truth (`benchmark/simulate.py` of flucut).
