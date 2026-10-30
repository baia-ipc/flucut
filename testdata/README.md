# Test data

The files of the program tests (`make test_programs`, see
`tests/run-case.sh`). The cases are listed in `cases.tsv`; the expected
digests of their results are in `expected/`.

## Reference

`ref-ha.fasta`: one real HA record, NCBI `NC_004908.1` (Influenza A
virus A/Hong Kong/1073/99 (H9N2), 1714 bases), downloaded with
`scripts/fetch_references.py`; its provenance (source URL and date) is
in `ref-ha.fasta.provenance.tsv`. The record lacks the conserved
termini: it starts with 11 of the 12 bases of the Uni12 terminus
(`GCAAAAGCAGG`) and ends after the coding sequence, without the Uni13
terminus.

The command was:

    python3 -I scripts/fetch_references.py --segment HA --subtype H5N1 \
        --max 1 --refseq-only --output ref-ha.fasta

It was asked for H5N1 and returned an H9N2 record; for these tests the
subtype does not matter, but the subtype filter of the script has to be
checked.

## Reads

`reads.fastq`: 12 reads **simulated** from the reference, not sequenced:

    python3 -I scripts/simulate-reads.py testdata/ref-ha.fasta \
        > testdata/reads.fastq

Each read is one situation (see the script): complete amplicons in both
orientations, flanking bases, substitutions, a read starting inside the
primer, amplicons internally primed by Uni12 and by Uni13 (at arbitrary
positions, not the site of a real artefact), a primer dimer, a read
without primers, a short amplicon and an amplicon with one primer only.

`barcode01/`: the same reads split into two files, the second one
gzipped, as a directory input (the output must be identical to that of
`reads.fastq`):

    head -n 24 testdata/reads.fastq > testdata/barcode01/part1.fastq
    tail -n +25 testdata/reads.fastq | gzip -n \
        > testdata/barcode01/part2.fastq.gz

`fastq_pass/`: a run-like directory, as written by MinKNOW, with two
barcode directories, each a group of its own (`barcode01`: the first 6
reads, `barcode02`: the other 6, gzipped), and a file beside them, which
is not read (a warning is output):

    mkdir -p testdata/fastq_pass/barcode01 testdata/fastq_pass/barcode02
    head -n 24 testdata/reads.fastq \
        > testdata/fastq_pass/barcode01/part1.fastq
    tail -n +25 testdata/reads.fastq | gzip -n \
        > testdata/fastq_pass/barcode02/part2.fastq.gz
    head -n 4 testdata/reads.fastq > testdata/fastq_pass/beside.fastq

`adapters.fastq`: the first 4 reads, each with the adapters of Oxford
Nanopore added (native barcode NB01 at the start, the reverse complement
of its bottom adapter at the end), once without its first 10 bases (the
read begins inside or after the primer: the adapter is what can be
found) and once with the ligation adapter instead at the start, no base
removed (see the script); for the cases `adapters` and `no_adapters`
(`--adapters none`):

    python3 -I scripts/add-adapters.py testdata/reads.fastq \
        > testdata/adapters.fastq
