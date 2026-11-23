<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)"
            srcset="docs/logo/flucut-banner-dark-640.png">
    <img alt="flucut" src="docs/logo/flucut-banner-640.png" width="480">
  </picture>
</p>

# flucut

Primer trimming for influenza A whole-genome sequencing.

flucut removes the primer-derived bases from the reads of the influenza A
multisegment RT-PCR with the universal Uni12/Uni13 primers, before the
assembly of the segments (e.g. with IRMA). It is not specific to a
sequencing platform: the error rates are parameters (the defaults are
those of Oxford Nanopore reads), and Illumina reads are handled as well.

## Why

The universal primers bind the conserved termini of the segments, which
are the same in all eight segments and in all subtypes. The bases of a
read which lie under a primer therefore come from the primer and not
from the template: a difference between the template and the primer at
such a position is overwritten, and a variant seen there is an artefact.

The Uni12 primer is short: used without a tail it is only the twelve
bases of the conserved terminus, with a degenerate position. A sequence
that similar also occurs inside the segments, where the primer can
anneal and start an amplicon in the middle of a segment. The reads of
such an amplicon begin with primer-derived bases which are not the
template there, and the primer writes its own sequence into the
consensus: a variant which looks real and is not, reproducible from run
to run, and present in the public databases. This is what flucut is
built around; see `docs/internal-priming.md`.

IRMA does not trim primers (it matches a single adapter sequence,
exactly or with one mismatch), so the trimming has to be done before it.

## What it does

- finds the Uni12/Uni13 primers at both ends of a read and in both
  orientations, whole or partial (a read can begin inside a primer),
  with errors (the defaults are those of Nanopore reads);
- removes, by position and not by difference, every primer-derived base
  (the default): this removes internal artefacts as well, as
  the primer bases of an internally primed read are at its end, like
  those of any other read. The consensus then lacks the 12 and 13 bases
  of the conserved termini, which would have been the primer's sequence,
  not measured sequence;
- keeps the conserved termini, for complete ends (the default): because
  flucut knows which primers were used, it knows where inside a primer
  the tail ends and the terminus begins. This mode uses references of the
  segments: the bases after the primer of each read are located in them,
  and a read which starts at an interior position of a segment (internally
  primed) loses its whole primer, the others only the tail. By default the
  references are a **built-in flu A set** (complete segments of several
  subtypes), downloaded from NCBI on first use and cached; `-r <file>`
  uses your own instead, and `-r off` turns references off (every
  primer-derived base removed). With several references of a segment
  (e.g. of several subtypes), each barcode is classified against the
  reference of each segment closest to its reads;
- with the references, in either mode, measures internal priming: which
  segments, which sites, what share of the reads (the references must be
  complete segments, with the UTRs, not coding sequences; flucut warns
  otherwise);
- a primer is accepted as the primer of an end only near that end of
  the read: a match deeper inside is not the primer which started the
  amplicon, and trimming there would remove template. Those matches are
  reported instead, as the sites at which the primer can anneal and as
  the junctions of chimeras;
- at a read end without a primer (a broken molecule), removes the
  sequencing adapters of Oxford Nanopore (native and rapid barcoding,
  ligation), and a primer remnant next to them;
- reads which are not amplicons are counted and reported: primer dimers,
  reads without a primer, reads with a primer at one end only.

## The primers

The conserved termini and the primers built into flucut are those of the
publications, quoted and not derived; each sequence is documented with
its source in `src/primers.h`. The built-in schemes are:

- `uni12-mbt13` (the default): the untailed degenerate Uni12 primer and
  the tailed Uni13 primer, the combination used in practice;
- `mbt`: the tailed primers of the multisegment RT-PCR;
- `uni`: the conserved termini used as primers, without a tail.

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

The main options (`flucut -h` lists all of them, with their defaults):

- `-r <references.fasta>`: classify each read end as terminal or
  internally primed against the segment references. The mode follows
  from this: **with** `-r`, the conserved termini are kept at the ends
  confidently classified as terminal (complete segment ends) and the
  whole primer is removed everywhere else (an internal or an
  unclassified end: cut when unsure, so a primer base never reaches the
  consensus unless the end is a known terminus); **without** `-r`, every
  primer-derived base is removed. The references are indexed by their
  k-mers, so that the time hardly depends on their number
  (`--exhaustive-search` aligns each read end to every reference, to
  check the index);
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
- `-p`: the platform: `nanopore` (default; reads of whole amplicons)
  or `illumina` (reads of fragments of the amplicons: the primer is at
  the read start; window 10, minimum length 1, no adapter search);
- `-w`: a primer must lie within this many bases of its end of the read
  (default 100; 10 with `-p illumina`);
- `-m`: reads shorter than this after trimming are not written
  (default 50; 1 with `-p illumina`, so that the two files of paired
  reads keep the same reads; primer dimers are never written);
- `-j`: worker threads (default 4), each trimming batches of 1000 reads
  in turn, written in the order of the input: the output does not
  depend on it;
- `--compress-level`: the gzip compression level of the outputs, 1 to
  12 (default 6; higher: smaller and slower; the outputs are compressed
  with libdeflate).

Outputs:

- `<name>.fastq.gz`: the trimmed reads, in their original orientation,
  with their titles; unless `--no-tags` is given, tab-separated SAM-like
  tags are appended to each title: `fs` (strand of the amplicon: `+`,
  `-`, or `.` if no primer was found), `fl`/`fr` (the primer found at the
  left/right end), `tl`/`tr` (bases removed at the start/end of the
  read), `ni` (primer matches inside the read), `sp` (part `k/n` of a
  read which was split, see below); with references also
  `ol`/`or` (origin of the end: terminal, internal or unknown),
  `ml`/`mr` (where the bases after the primer mapped:
  reference:position:strand) and `ct` (the category of the read: `full`,
  `internal`, `deletion` or `other`, see below). IRMA accepts the tags (same consensus
  with and without them, see `benchmark/RESULTS.md`); `--no-tags` writes
  the titles unchanged. A read which contains several amplicons joined
  end to end is split, and each part is written as its own record, named
  `<read>_1`, `<read>_2`, ...; `--no-split` turns this off.
- with `--by-category`, four files per sample instead of
  `<name>.fastq.gz`, one for each category of reads (the colours of the
  coverage plot of the report):

  | file | reads |
  |---|---|
  | `<name>.full.fastq.gz` | full templates: the amplicon started at a segment end |
  | `<name>.internal.fastq.gz` | internal priming: it started inside the segment |
  | `<name>.deletion.fastq.gz` | internal deletions: the two ends of a segment joined |
  | `<name>.other.fastq.gz` | no 5' primer, or not located on the references |

  ```
  flucut --by-category -o trimmed barcode01/
  ```

  The four files together hold the same reads as `<name>.fastq.gz`
  would. The category needs the references (not `-r off`).
- `summary.tsv`: one line per output, with the reads and bases read and
  written, the reads by number of primers found, those on the minus
  strand, the dimers, the reads too short, the reads with internal
  primer matches, and, with references, the left and right ends by
  origin, the complete amplicons shorter than 90% of their reference
  (`short_complete`) and the reads written by category (`out_full`,
  `out_internal`, `out_deletion`, `out_other`); `NA` without references.
- `lengths.tsv`: the lengths of the complete amplicons (the template
  between the two primers), in bins of 10 bases: output, group (the
  reference to which the read end mapped, `not_located`, or `all`
  without references), start of the bin, reads. It is a check
  independent of the classification of each read end: internal priming
  (and deletions, e.g. in defective interfering RNAs) shows as complete
  amplicons much shorter than their segment. The log names, with
  references, each reference with such amplicons and the length where
  most of them are; without references, it lists the peaks of the
  lengths.
- `flucut_report.html`: the HTML report, one self-contained file
  (inline SVG plots, light and dark themes), made to be read at a
  glance:
  - a **run overview**: the samples by state, and a table with one row
    per sample: complete segments (e.g. `7/8`), subtype, reads, N50 and
    quality, a cell for each segment (PB2 to NS) with the share of it at
    a depth of 50 or more, the segments with most internal priming and
    most internal deletions, and the reads by origin; and the median
    depth of each segment in each sample, as a plot;
  - the **HA in all samples**, side by side;
  - for each **sample**: its state and what deserves attention, the
    subtype, the reads by primers found and by origin, the internal
    priming sites, the bases by segment, a strip for each segment
    (numbers, coverage by read origin where internal priming shows as a
    step up, lengths of the complete amplicons), and the lengths and
    qualities of the reads before and after the trimming, with the
    flow cell, run and basecalling model found in the read titles.

  See the Reports chapter of the manual for how to read it.
- `flucut_dashboard.html`: the same content, one sample on the screen
  at a time, to follow the run (below).
- `flucut_summary.json`: the same content for programs, in a versioned
  format (`docs/json_summary.md`).

A run with several samples can be followed while it is in progress:
flucut writes the dashboard before the first sample and after each one,
and says where it is (`Dashboard: out/flucut_dashboard.html (open it in
a web browser to follow the run)`). The page reloads itself while the
run is in progress; at the end it stays as the dashboard of the
finished run. The report is a separate page with all the samples one
after the other, for printing or archiving. In the dashboard the strips
of the segments are folded (a click opens one), and the HA tab starts
with a table of the HA of every sample. The messages show which
sample is being trimmed (`Trimming run1/barcode01.fastq.gz (1 of 3)...`).
If such a run is interrupted, running it again with `--resume` (same
inputs, options and output directory) skips the samples already
trimmed and trims only the others; on a complete run, `--resume` trims
nothing and only writes the tables and reports again. The state of each
sample is kept in `outdir/.flucut-state/`.

## With IRMA

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

## Checking the predictions

The report of flucut is built from probes, not from mappings. To see how
close it is for a run, `posthoc/flucut-posthoc` maps every read to the
segment references and measures the same things from the mappings (the
origin of each amplicon end, the internal priming sites, the amplicon
lengths, the coverage by read origin), next to what flucut predicted and
to what IRMA assembled. It is a separate program, in Python; see
`posthoc/README.md`.

## Live mode

`flucut --live -o out /data/run1/fastq_pass` follows a MinKNOW run while
it is sequencing: the files are
trimmed as they are written (once not modified for `--min-age`, default
60 s, and complete), appended to the outputs of their samples, and the
tables, the dashboard and the JSON summary are updated after each file
(the dashboard reloads itself). It stops at the end of the run (the `final_summary_*.txt` of
MinKNOW), after `--timeout` without new files, or with Ctrl-C; with
`--resume`, an interrupted live run continues where it stopped.

## Manual

The user manual is `manual/flucut-manual.pdf` (sources in
`docs/manual`; `make manual` rebuilds it, with Sphinx and LaTeX).

The methods are described algorithm by algorithm, as pseudocode with
explanations and figures, in `docs/methods/methods.md`, also as
`manual/flucut-methods.pdf` (`make methods` rebuilds it; this needs
Sphinx with `myst-parser`, matplotlib and XeLaTeX).

## Building

Only a C compiler and `make` are needed: the libraries used (libcbaia,
edlib, libyaml, zlib, libdeflate) are included in `ext` (see
`ext/README.md`).

```
make
make test
```

`make SYSTEM_LIBS=1` uses the zlib and libyaml of the system instead of
the included ones (libdeflate, which compresses the output, is always
the included one).

## State

Early development. Implemented and tested: the built-in primers and
schemes, the search of the primers (degenerate bases, partial primers
at the read ends), the trimming of a read, the classification of the
read ends against references, and the command line program with its
gzipped FASTQ output, summary table, HTML report and JSON summary, the
groups of input files (MinKNOW barcode directories), the worker
threads, the adapter search, `-e auto` and the live mode.

## License

ISC (see `LICENSE`).
