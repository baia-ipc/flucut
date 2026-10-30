# The tools compared

Each tool is configured as well as its documentation allows; the exact
command of every run is recorded here and in `RESULTS.md`, together with
the version used. Where a tool is not meant for a task, this is stated
instead of running it badly: a comparison is only useful if each tool is
given its best configuration.

## Barbell

Rust; pattern-aware demultiplexing and trimming of Nanopore reads.
Beeloo R, Groot Koerkamp R, Jia X, Broekhuizen-Stins MJ, van Ijken L,
Dutilh BE, Zomer AL. *Barbell resolves demultiplexing and trimming
issues in Nanopore data.* Bioinformatics (2026),
https://doi.org/10.1093/bioinformatics/btag349 ;
https://github.com/rickbeeloo/barbell

The reference to beat. It reports far fewer trimming errors than
`dorado trim` and Flexiplex, and it is fast. Its model is close to the
one needed here: a query FASTA holds the full expected sequences, of
which it takes the shared prefix as the adapter flank and the shared
suffix as the primer flank, the differing middle being the barcode.

For the multisegment RT-PCR this gives, without barcodes, one entry per
end, and with barcoded primers one entry per sample:

    >sample
    <adapter><barcode><primer>

    barbell annotate -q left.fasta,right.fasta -b Ftag,Rtag \
      -i reads.fastq -o anno.tsv -t <threads>
    barbell inspect -i anno.tsv
    barbell filter -i anno.tsv -f filters.txt -o filtered.tsv
    barbell trim -i reads.fastq -a filtered.tsv -o out

with a pattern such as

    Ftag[fw, @left(0..250), >>]__Rtag[<<, rc, @right(0..250)]

Points to check in the comparison, honestly, in both directions:

- Barbell trims the whole matched query, so it removes the conserved
  terminus together with the tail. That is what flucut does by default
  too, so on total primer removal the two should be close; the
  difference to measure is in the **boundary**: how many template bases
  are removed with the primer, and how many primer bases are left.
- A read which begins *inside* a primer (very common at the ends of
  amplicons) has no complete flank. Whether Barbell recovers these, and
  how many template bases it then removes, is one of the main things to
  measure.
- Barbell is a demultiplexer first: with a single sample and no barcode,
  the flank extraction degenerates (prefix and suffix of a single
  sequence). The fair configuration for that case must be found and
  documented, not assumed.
- Its `--use-extended` mode searches for fusion points and other
  artefacts (about 3x slower) and is recommended when consensus quality
  matters, which is exactly our case: it must be included, as well as
  the default mode.
- Barbell has its own simulation and comparison harness
  (`benchmarks/src/simulations`, `benchmarks/src/compare`). Where our
  metrics can be made to match theirs, they are, so that the numbers can
  be compared with the published ones.

## IRMA-core trimmer

The trimmer used today in the CDC MIRA-NF pipeline
(https://github.com/CDCgov/MIRA-NF), with `custom_primers`,
`primer_kmer_len` and `primer_restrict_window`. The closest tool in
purpose: it is meant for exactly this step, before IRMA.

## dorado trim

Oxford Nanopore's basecaller trimming, with `--primer-sequences` for
custom primers.

## Porechop / Porechop_ABI

Adapter trimming for Nanopore reads, with the primers added as custom
adapters.

## cutadapt

General-purpose adapter trimming; linked adapters and anchored forms are
used where they fit, with the error rate set for Nanopore reads.

## BBDuk

k-mer based trimming (`ktrim`, `k`, `hdist`, `restrictleft`,
`restrictright`).

## fastp

Mainly for Illumina (the second phase); its adapter options are used for
the paired-end data.

## Alignment-based trimming

`ivar trim` and `samtools ampliconclip` trim by primer position on a
mapped BAM, given a primer BED. They need a mapping step, which the
others do not, and they are the natural comparison for "trim by
position" rather than by sequence. They are run on reads mapped to the
segment references, with the primer intervals of the conserved termini.

## IRMA itself

IRMA has no primer trimming: it has a single `ADAPTER` sequence, matched
exactly or with one mismatch. It is therefore not a tool in the
comparison, but the consumer of the output: the effect of each tool on
the IRMA consensus, at the segment ends, is one of the measures.
