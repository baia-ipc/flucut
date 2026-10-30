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
  (the default): this removes the internal priming artefact as well, as
  the primer bases of an internally primed read are at its end, like
  those of any other read. The consensus then lacks the 12 and 13 bases
  of the conserved termini, which would have been the primer's sequence,
  not measured sequence;
- optionally keeps the conserved termini, for complete ends: because
  flucut knows which primers were used, it knows where inside a primer
  the tail ends and the terminus begins. This mode requires the
  references of the segments (`scripts/fetch_references.py` downloads
  them from NCBI): the bases after the primer of each read are located
  in them, and a read which starts at an interior position of a segment
  (internally primed) loses its whole primer, the others only the tail;
- with the references, in either mode, measures internal priming: which
  segments, which sites, what share of the reads;
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

- `-T primer` (default): remove every primer-derived base;
  `-T terminus`: keep the conserved termini, which requires
  `-r <references.fasta>`;
- `-r <references.fasta>`: classify each read end as terminal or
  internally primed against the segment references; with `-T primer`
  the classification is only reported; the references are indexed by
  their k-mers, so that the time hardly depends on their number
  (`--exhaustive-search` aligns each read end to every reference, to
  check the index);
- `-e`: errors allowed in a primer (a fraction of its length, a number
  of bases, or `auto`, the default: 0.2, or 0.3 if the primers of the
  first 1000 reads have at least 4% errors);
- `-w`: a primer must lie within this many bases of its end of the read
  (default 100);
- `-m`: reads shorter than this after trimming are not written
  (default 1; primer dimers are never written);
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
  read), `ni` (primer matches inside the read); with references also
  `ol`/`or` (origin of the end: terminal, internal or unknown) and
  `ml`/`mr` (where the bases after the primer mapped:
  reference:position:strand). IRMA accepts the tags (same consensus
  with and without them, see `benchmark/RESULTS.md`); `--no-tags` writes
  the titles unchanged.
- `summary.tsv`: one line per output, with the reads and bases read and
  written, the reads by number of primers found, those on the minus
  strand, the dimers, the reads too short, the reads with internal
  primer matches, and, with references, the left and right ends by
  origin and the complete amplicons shorter than 90% of their reference
  (`short_complete`; `NA` without references).
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
- `flucut_report.html`: the HTML report, self-contained (inline SVG
  plots, light and dark themes): the run parameters and, for all reads
  and for each group, the reads at each step, the primers found at each
  end, with references the internal priming sites with their reads, the
  complete amplicons shorter than 90% of their reference and the
  histograms of the lengths.
- `flucut_summary.json`: the same content for programs, in a versioned
  format (`docs/json_summary.md`).

## Live mode

`flucut --live -o out /data/run1/fastq_pass` follows a MinKNOW run while
it is sequencing (as the real-time mode of Lumipore): the files are
trimmed as they are written (once not modified for `--min-age`, default
60 s, and complete), appended to the outputs of their samples, and the
tables and reports are updated after each file (the HTML report reloads
itself). It stops at the end of the run (the `final_summary_*.txt` of
MinKNOW), after `--timeout` without new files, or with Ctrl-C; with
`--resume`, an interrupted live run continues where it stopped.

## Manual

The user manual is `manual/flucut-manual.pdf` (sources in
`docs/manual`; `make manual` rebuilds it, with Sphinx and LaTeX).

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
