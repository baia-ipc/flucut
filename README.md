<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)"
            srcset="docs/logo/flucut-banner-dark-640.png">
    <img alt="flucut" src="docs/logo/flucut-banner-640.png" width="480">
  </picture>
</p>

# flucut

flucut is a comprehensive trimming, filtering and QC-tool for influenza A and B whole-genome sequencing. 

In a very fast pass over the reads, flucut will:
- classify the amplicons, including:
  - those from the entire template, which are the intended ones
  - those originating from inexact internal priming, which can cause wrong base calls in
    assembly programs such as IRMA
  - those containing internal deletions, e.g. originating from DIPs, which in some case
    cause misassemblies over the junction
  - those from fragmented amplicons
  - reads containing the concatenation of multiple amplicons
- remove adapters and other library-derived spurious sequences beyond the primers
- remove the primer-derived bases in internal priming segments 
- split reads containing multiple amplicons joined end-to-end
- estimate the coverage plots for each segment accurately before any read mapping is done
- classify the subtype of influenza A
- report reads of a second HA or NA subtype in a sample (a mixed infection, a contamination
  or barcode crosstalk: the reads do not say which)
- analyze presence of patterns typical of low/high pathogenicity
- read the amino-acid markers of host adaptation, receptor binding and antiviral
  resistance (PB2 627 and 701, PA 38, NA 275, HA 226/228, M2 31, the NA stalk,
  the NS1 C-terminus and PB1-F2) straight from the reads, with the share of the
  reads carrying each residue, so that a minority a consensus would hide is seen

It is mainly targeted to Nanopore reads, but also supports most functionality on Illumina reads as well.

The classification of the amplicon type, subtype and markers requires reference sequences.
These can be (1) a built-in set, which is downloaded from NCBI on first use (2) the seed
sequences used by IRMA (3) a user-provided reference dataset.

## Building

Only a C compiler and `make` are needed
(the source code of all libraries used, i.e. libcbaia,
edlib, libyaml, zlib, libdeflate is included in the `ext` directory).

To compile just run:
```
make
```

The executable is contained in the ``bin`` directory and does not require any other file,
thus it can be moved or linked to any directory the user desires (e.g. one specified
in the PATH environmental variable).

## Usage

```
flucut -o trimmed sample1.fastq.gz barcode02/
```

Each input is a FASTQ file (optionally gzipped), an unaligned BAM file,
or a directory of such files (e.g. a barcode directory of a MinKNOW
run, whose files are merged); a directory with `barcode*`
subdirectories (the `fastq_pass` directory of a MinKNOW run) gives one
group of files for each of them (the files beside them are not read).
The trimmed reads of each group are written to
`<outdir>/<name>.fastq.gz`, one file per sample, ready for IRMA.
`<name>` is the file name without its suffix, or the directory name. Two
groups with the same name are an error.

`flucut -h` lists all options, with their defaults.

Main options:
- `-p`: the platform: `nanopore` (default; reads of whole amplicons)
  or `illumina` (reads of fragments of the amplicons: the primer is at
  the read start; window 10, minimum length 1, no adapter search);
- `-j`: worker threads (default 4), each trimming batches of 1000 reads
  in turn, written in the order of the input: the output does not
  depend on it;

Other options:
- `-r <references.fasta>`: classify each read end as terminal or
  internally primed against the segment references. The mode follows
  from this: **with** `-r`, the conserved termini are kept at the ends
  confidently classified as terminal (complete segment ends) and the
  whole primer is removed everywhere else (an internal or an
  unclassified end: cut when unsure, so a primer base never reaches the
  consensus unless the end is a known terminus); **without** `-r`, every
  primer-derived base is removed. The references are indexed by their
  k-mers, so that the time hardly depends on their number
- `--virus a|b|both` says which types a run may hold. The default is `both`,
   which incurs in a small performance loss, of a fraction of a second per run.

- `-s <scheme>`: the primer scheme, i.e. the primers searched at the read
  ends (see "Supported primers" below). The default, `auto`, chooses the
  scheme of influenza A from the primer tails found in the first reads; where
  they do not single one out, the scheme of `--virus` is used: `uni12-mbt13`
  for `a`, `ibv-ga2` for `b` and `uni12-mbt13+ibv-ga2` for `both`.

References from the seeds of IRMA:
- the seed sequences of IRMA hold the coding sequence of each segment only,
  without the untranslated regions where the primers bind, so they cannot be
  given to `-r` as they are: the genuine segment ends would be classified
  as internal priming
- `scripts/irma_references.py` builds a reference set from them: to each
  seed it attaches the untranslated regions of a complete record of the same
  segment (for HA and NA, of the same subtype), taken from the built-in set
  ```
  python3 scripts/irma_references.py --irma /opt/flu-amd/IRMA_RES --save-config
  ```
- `--save-config` stores the paths of IRMA and of the set written in the
  configuration file of flucut (`$XDG_CONFIG_HOME/flucut/config`, else
  `~/.config/flucut/config`): a run without `-r` then uses this set instead
  of the built-in one (`-r <references.fasta>` and `-r off` still override it)
- `--from-config` builds the set again from the stored paths, e.g. after
  IRMA updates its seeds
- this set is meant for trimming only, never as a reference for an assembly:
  the coding sequence and the untranslated regions of each sequence come from
  two different viruses

Outputs:
- `<name>.fastq.gz`: the trimmed reads, in their original orientation,
  with their header lines and SAM tags appended to them (see below)
- `flucut_dashboard.html`: a self-updating dashboard page, which can be
followed during the program run, giving program run diagnostics and ETA,
metadata of the sequences and all the QC information collected, e.g.
 completeness of the segments, subtype, number of reads, N50, quality, length histogram,
 median depth, classification in amplicon categories, e.g. internal priming and internal deletions, coverage plots; 
it includes an overview of the sequencing, a page for each sample
and a page focused on the HA segments and their cleavage site
- `flucut_report.html`: the same content in a printable form,
constructed at program completion
- `flucut_summary.json`: the same content in JSON for programmatic access

## Interplay with MinKNOW (Live mode)

Live mode:
- `--live` follows a MinKNOW run while it is still sequencing, so no time is lost:
the files are processed as they are completely written
- `--min-age`: time after which a sequencing file is considered final
(completely written) and thus processed by flucut (default: 60 s)
- `--timeout`: time without new files appearing after which flucut
live automatically stops (in completed runs, once the `final_summary_*.txt`
file written by MinKNOW appears, flucut is stopped)
- `--resume`: resume an interrupted live run, where it stopped

## Interplay with IRMA

flucut trims the primers; IRMA assembles the segments. IRMA does not trim
the universal primers (it matches a single adapter, exactly or with one
mismatch), so without flucut the primer-derived bases, including those of
the internally primed reads, are written into the consensus as false
variants (e.g. the HA artefact at the Uni12-like site). Run flucut first,
then IRMA on the trimmed reads, one sample at a time.

Nanopore (the `FLU-minion` module of IRMA):

```
flucut -o trimmed /data/run1/fastq_pass
for fq in trimmed/*.fastq.gz; do
  sample=$(basename "$fq" .fastq.gz)
  IRMA FLU-minion "$fq" "irma/$sample"
done
```

`flucut` writes one trimmed `trimmed/<barcode>.fastq.gz` per sample (plus
the tables and the report); IRMA reads one of them at a time and writes
the per-segment consensus to `irma/<barcode>/` (`A_HA_*.fasta`,
`A_NA_*.fasta`, ...). IRMA accepts the tags flucut appends to the read
titles (same consensus with and without them; `--no-tags` removes them).

Illumina: trim with `flucut -p illumina` and assemble with the paired
`FLU` module of IRMA instead of `FLU-minion`.

## Output by category

- with `--by-category`, six files per sample are written instead of
  `<name>.fastq.gz`, one for each category of reads
- here is the list of files:
  | file | reads |
  |---|---|
  | `<name>.full.fastq.gz` | full templates: the amplicon started at a segment end |
  | `<name>.internal.fastq.gz` | internal priming: it started inside the segment |
  | `<name>.deletion.fastq.gz` | internal deletions: the two ends of a segment joined |
  | `<name>.fragment.fastq.gz` | fragments: placed on a segment, but no primer at the 5' end (broken molecules, truncated reads) |
  | `<name>.unplaced.fastq.gz` | complete amplicons the references place on no segment |
  | `<name>.unknown.fastq.gz` | the other reads the references place on no segment: mostly not influenza |
- the first four are the colours of the coverage plot of the report: a
  read is written there only if the references say which segment it is.
  A read with a primer at both ends whose probes did not map goes to
  `unplaced`, even where its length is that of a segment, so that
  `<name>.full.fastq.gz` holds influenza reads only
- command line example
  ```
  flucut --by-category -o trimmed barcode01/
  ```
- the six files together hold the same reads as `<name>.fastq.gz`
  would.
- this option needs the references, i.e. it is incompatible with `-r off`

## Reads with an internal deletion, chimeras, reads on no segment

- a read with an internal deletion (e.g. of a defective interfering
  genome) is two flanks of a segment joined at a junction. By default it
  is written as two records, its two flanks, named `<read>_j1` and
  `<read>_j2` (tag `jf:Z:1/2`, `jf:Z:2/2`): a mapper then gets two
  reads which each lie in one place, instead of one read with a gap, and
  the reads of the deletion genomes no longer bias the consensus at the
  junctions. A flank is also cut at each short deletion (6 bases or
  more) in its 250 bases next to the junction, which many deletion
  genomes carry and which would otherwise enter the consensus: the read
  is then written as more pieces (`_j1`, `_j2`, `_j3`, ...). A read whose
  junction could not be located is written whole
- `--deletions keep` writes these reads whole, `--deletions filter`
  leaves them out (an assembly of the reads without the deletion
  genomes):
  ```
  flucut --deletions filter -o trimmed barcode01/
  ```
- a read which contains several amplicons joined end to end (a chimera)
  is cut into its amplicons by default; `--chimeras keep` writes it
  whole, `--chimeras filter` leaves it out
- `--unplaced filter` leaves out the complete amplicons the references
  place on no segment, `--unknown filter` the other reads they place on
  no segment (mostly not influenza)
- a read left out is still analyzed and counted: the reports describe
  all the reads of the sample, and give the reads left out (`filtered`)
- `cut` and `filter` need the references; with `-r off` the reads with
  an internal deletion are written whole

## Markers

Besides the HA cleavage site, flucut reads a few positions of the influenza A
proteins whose residue is worth knowing before anything is assembled, straight
from the reads: PB2 627 and 701 (replication in a mammalian host), PA 38
(baloxavir), NA 275 (oseltamivir), HA 226 and 228 (which receptors the virus
binds), M2 31 (the adamantanes), the length of the NA stalk, the C-terminus of
NS1 and the reading frame of PB1-F2.

Each position is found by a short conserved peptide directly beside it, which
is searched in the read itself, so no position in the sequence is used and a
virus whose numbering is shifted still gives the right residue. What is
reported is a distribution, not a call: the residue, the reads it was read
from, and the share of those reads carrying each value, so that a residue
present in a part of the sample is visible instead of being rounded away by a
consensus.

They are in the `markers` array of `flucut_summary.json`, in a column per
marker in `summary.tsv`, and in the **Markers** card of each sample of the
report and the dashboard. Clades and genotypes are deliberately not reported:
they are defined by nomenclatures which are revised, which a table compiled
into a program cannot follow. See `docs/markers.md` and the manual chapter.

## Evaluation suite

An evaluation suite named `flucut-posthoc` was written in Python, with the main
purpose of evaluating how valid are the heuristics of flucut, in comparison
with exact information collected through read mapping and through IRMA.

## Manual

A comprehensive user manual is provided as PDF in `manual/flucut-manual.pdf`
(sources in `docs/manual`; `make manual` rebuilds it, requiring Sphinx and LaTeX).

Additionally, a technical manual describes the algorithms employed in flucut
as pseudocode, with further explanations and figures, and is provided as PDF
in `manual/flucut-methods.pdf` (sources in `docs/methods`;
`make methods` rebuilds it, requiring Sphinx with `myst-parser`, matplotlib and XeLaTeX).

## What flucut knows about influenza

The sequences and tables flucut carries are not in the C sources. They are in
`metadata/`, one readable file per kind, which the build turns into the tables
of the program: `primers.yaml` (the primers, the conserved termini and tails,
the segments each binds, the source of each sequence, and the schemes),
`markers.yaml` (the markers and the anchor of the HA cleavage site),
`adapters.yaml` (the Nanopore kit adapters), `references.yaml` (the accessions
of the built-in reference sets), `segments.yaml` (the segments and the words
which name their genes) and `examples.yaml` (the real sequences the unit tests
use).

Adding a primer, a marker or a kit is a record in one of those files; no C
source changes with it. The generators are awk, so building still needs only a
C compiler and `make`. See the manual, *Installation*.

## Influenza A and influenza B

The trimming, the classification of the reads, the output by category,
the reports and the HA cleavage site work the same for both types. These
are of influenza A only: the subtype (a sample of influenza B is `B`, or
`B/Victoria` and `B/Yamagata` where the reference names the lineage), the
reads of a second HA or NA subtype, the markers, the choice of the
primer scheme from the tails (`-s auto` chooses among the schemes of
influenza A), and `flucut-posthoc`. Influenza B has been tested on
simulated reads of three genomes and on four public Illumina runs of
B/Victoria; no Nanopore run of influenza B yet. See the manual, *The
primers*.

## Supported primers

These are the supported primer schemes:
- `uni12-mbt13+ibv-ga2` (used when `-s auto` singles out no scheme): the
  primers of both types together, for a run which may hold samples of
  either;
- `uni12-mbt13`: the untailed degenerate Uni12 primer and the tailed
  Uni13 primer, the combination used in practice for influenza A;
- `mbt`: the tailed primers of the multisegment RT-PCR of influenza A;
- `uni`: the conserved termini of influenza A used as primers, without
  a tail;
- `opti`: the Opti1 primers of influenza A (Mena et al. 2016);
- `uni-inf`: the tailed universal primers Uni-12/Inf and Uni-13/Inf1 of
  influenza A;
- `ibv-ga2`: the universal IBV-GA2 primer cocktail of influenza B
  (Zhou et al. 2014), six pairs for the eight segments.

Several schemes joined by `+` are one scheme with the primers of each
(e.g. `opti+ibv-ga2`). `-s auto` (the default) chooses the scheme of
influenza A from the tails found directly in front of the conserved
termini in the first reads; the manual (*The primers*) says how, and
what happens when no single scheme fits.

## SAM Tags

Unless `--no-tags` is given, tab-separated SAM-like
tags are appended to each read header line: `fs` (strand of the amplicon: `+`,
`-`, or `.` if no primer was found), `fl`/`fr` (the primer found at the
left/right end), `tl`/`tr` (bases removed at the start/end of the
read), `ni` (primer matches inside the read), `sp` (part `k/n` of a
read which was split, see below); with references also
`ol`/`or` (origin of the end: terminal, internal or unknown),
`ml`/`mr` (where the bases after the primer mapped:
reference:position:strand) and `ct` (the category of the read: `full`,
`internal`, `deletion`, `fragment`, `unplaced` or `unknown`, see below).

## Advanced user and developer options

Advanced user options:
- `--segment <s>`: for runs of one segment alone (e.g. `--segment HA`):
  classify the read ends against the references of that segment only
  (faster, no chance matches in the other segments);
- `--report-segment <s>`: show only one segment in the HTML report, by
  number or gene (e.g. `--report-segment HA`); the trimming does not
  change;
- `-e`: errors allowed in a primer: a percentage of its length (`-e 20%`,
  or as a fraction, `-e 0.2`), a number of bases (`-e 3`), or `auto`,
  the default: 20%, or 30% if the primers of the first 1000 reads have
  at least 4% errors;
- `-w`: a primer must lie within this many bases of its end of the read
  (default 100; 10 with `-p illumina`);
- `-m`: reads shorter than this after trimming are not written
  (default 50; 1 with `-p illumina`, so that the two files of paired
  reads keep the same reads; primer dimers are never written);
- `--compress-level`: the gzip compression level of the outputs, 1 to
  12 (default 6; higher: smaller and slower; the outputs are compressed
  with libdeflate).
- `--resume`: resumes an interrupted run, or reconstruct the report

Developer option:
- `--exhaustive-search` aligns each read end to every reference, to
  check the index
- `--no-tags`: do not add SAM tags to the read header lines

Advanced build settings:
- `make SYSTEM_LIBS=1` uses the zlib and libyaml of the system instead of
the included ones

## License

ISC (see `LICENSE`).
