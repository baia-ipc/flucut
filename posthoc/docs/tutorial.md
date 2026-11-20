# Tutorial

## One sample, before and after flucut

You have sequenced a sample, trimmed it with flucut and assembled it with
IRMA, once from the untrimmed reads and once from the trimmed ones:

    run1/barcode01.fastq.gz               the reads of the sequencer
    trimmed/barcode01.fastq.gz            written by flucut
    trimmed/flucut_summary.json
    trimmed/flucut_report.html
    irma/barcode01_IRMA/                  IRMA on the untrimmed reads
    irma/barcode01_flucut_IRMA/           IRMA on the trimmed reads

Run:

    ./flucut-posthoc -v -j 16 -o posthoc/barcode01 \
        -i irma/barcode01_IRMA,irma/barcode01_flucut_IRMA \
        -f trimmed/barcode01.fastq.gz \
        -s trimmed/flucut_summary.json -R trimmed/flucut_report.html \
        run1/barcode01.fastq.gz

`-v` reports the progress. The group of the sample in the outputs of
flucut is found from the name of the reads file (`barcode01`); if they
differ, name the group with `-g`, and the sample with `-n`.

Open `posthoc/barcode01/posthoc_report.html` and read it from the top:

1. **The tiles and "What the reads are"**: how many reads are influenza
   at all. A run can be mostly amplicons of something else (the primers
   also amplify host sequences); these are grouped under "Amplicons
   which are not influenza", with a sequence to search in a database.
2. **References**: which reference each segment was compared with, and
   whether the consensus of IRMA covers it. A consensus much shorter
   than the reference was built from genomes with a deletion; one with
   primers inside was built from amplicons joined end to end.
3. **Internal priming sites**: where the primers annealed inside the
   segments, and how many reads flucut reports there.
4. **Coverage by read origin**: the measured depth at the left, the plot
   of flucut at the right, on the same scale.
5. **flucut compared with the mapping**: how often flucut's call for a
   read end is confirmed, and whether its trimming removed the primer.

## The reads of flucut only

Without the untrimmed reads, give the reads written by flucut as the
reads:

    ./flucut-posthoc -o posthoc/barcode01 -i irma/barcode01_flucut_IRMA \
        -s trimmed/flucut_summary.json -R trimmed/flucut_report.html \
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
          -s trimmed/flucut_summary.json -R trimmed/flucut_report.html \
          $fq
    done
    ./flucut-posthoc-overview -o posthoc posthoc/*/posthoc_summary.json

`posthoc/overview.html` has one line per sample, with the main measures
and the agreement of flucut, and links to the reports;
`posthoc/overview.tsv` is the same table for programs.

To compare one segment across the samples, add `--segment`:

    ./flucut-posthoc-overview --segment HA -o posthoc posthoc/*/posthoc_summary.json

`posthoc/coverage_HA.html` then shows, for each sample, the coverage of
HA measured from the mapping next to the one flucut predicts, on the
same scale, and a table of how well the two agree.

## Your own references

By default the references are the set which flucut downloads
(`~/.cache/flucut/flu-a-references-v1.fasta`); run flucut once to get
it, or give a FASTA file of complete segments with `-r`. As in flucut,
the segment of a record is read from its description (`segment 4`); from
several records of a segment the closest to the sample is used. Records
without a segment number are all used.
