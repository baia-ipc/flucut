# Tutorial

## One sample, before and after flucut

You have sequenced a sample, trimmed it with flucut and assembled it with
IRMA, once from the untrimmed reads and once from the trimmed ones:

    run1/barcode01.fastq.gz               the reads of the sequencer
    trimmed/barcode01.fastq.gz            written by flucut
    trimmed/flucut_summary.json
    irma/barcode01_IRMA/                  IRMA on the untrimmed reads
    irma/barcode01_flucut_IRMA/           IRMA on the trimmed reads

Run:

    ./flucut-posthoc -v -j 16 -o posthoc/barcode01 \
        -i irma/barcode01_IRMA,irma/barcode01_flucut_IRMA \
        -f trimmed/barcode01.fastq.gz \
        -s trimmed/flucut_summary.json \
        run1/barcode01.fastq.gz

`-v` reports the progress. The group of the sample in the outputs of
flucut is found from the name of the reads file (`barcode01`); if they
differ, name the group with `-g`, and the sample with `-n`.

Open `posthoc/barcode01/posthoc_report.html` and read it from the top:

1. **The two verdicts**: whether the subtype flucut assigned is the one
   the mapping shows, and how many checks are outside their limits (the
   chips beside them name those checks).
2. **Subtype**: for HA and NA, the type and the reference of the
   mapping, of flucut and of IRMA, with the reads behind the mapping.
3. **Agreement by measure, by segment**: every check with its value,
   what it is a share of, and its limits. `!` and `✕` are where to look.
4. **The sample, measured**: how many reads are influenza at all (a run
   can be mostly amplicons of something else), what the amplicons are,
   and the internal priming sites with the reads flucut reports there.
5. **Segments**: one strip per segment, with the measured coverage
   above the one flucut predicts, on the same scale.
6. **All numbers**: open a table when you need the exact counts.

## The reads of flucut only

Without the untrimmed reads, give the reads written by flucut as the
reads:

    ./flucut-posthoc -o posthoc/barcode01 -i irma/barcode01_flucut_IRMA \
        -s trimmed/flucut_summary.json \
        trimmed/barcode01.fastq.gz

Their tags are compared with their mapping as above. The primers are no
longer in the reads, so an end has a primer only where flucut kept the
conserved terminus, and the trimming cannot be checked.

## A whole run

Run each sample into its own directory, then make the overview:

    for fq in run1/*.fastq.gz; do
      name=$(basename $fq .fastq.gz)
      ./flucut-posthoc -j 32 -o posthoc/$name \
          -i irma/${name}_IRMA,irma/${name}_flucut_IRMA \
          -f trimmed/$name.fastq.gz \
          -s trimmed/flucut_summary.json \
          $fq
    done
    ./flucut-posthoc-overview -o posthoc posthoc/*/posthoc_summary.json

`posthoc/overview.html` has one line per sample: whether flucut agrees
with the mapping, the subtype of flucut, of the mapping and of IRMA,
the main measures and the checks, coloured by their limits, and a link
to the report; `posthoc/overview.tsv` has the same values for programs.

To compare one segment across the samples, add `--segment`:

    ./flucut-posthoc-overview --segment HA -o posthoc posthoc/*/posthoc_summary.json

`posthoc/coverage_HA.html` then shows the strip of HA of each sample
(numbers, coverage measured and predicted, amplicon lengths) and a
table of how well the two depths agree.

## Your own references

By default the references are the set which flucut downloads
(`~/.cache/flucut/flu-a-references-v1.fasta`); run flucut once to get
it, or give a FASTA file of complete segments with `-r`. As in flucut,
the segment of a record is read from its description (`segment 4`); from
several records of a segment the closest to the sample is used. Records
without a segment number are all used.
