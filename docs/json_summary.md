# JSON summary

flucut always writes a machine-readable summary of the run to
`outdir/flucut_summary.json`, with the same information as the HTML
report (`outdir/flucut_report.html`). The file is written to a temporary
file and then renamed, so that it can be read at any time. The format is
that of the summary of Lumipore, with the content of flucut.

## Format (version 1)

Top level:

| key | type | description |
|---|---|---|
| `format` | integer | version of this format (1) |
| `program`, `version` | string | `"flucut"` and its version |
| `updated` | string | time of writing (ISO 8601) |
| `running` | boolean | a run in progress, still followed (always `false` at the moment) |
| `command` | string | the command line |
| `scheme` | object | the primer scheme: `name`, and the names of the primers of the `left` (Uni12) and `right` (Uni13) ends, arrays of strings, in the order of the scheme |
| `trim_policy` | string | `"primer"` or `"terminus"` (`-T`) |
| `references` | object or null | with `-r`: `file` (the path given), `sequences` (array of `name`, `length`, `segment` (1-8, from the FASTA description, see the manual) and `gene` (e.g. `HA`), `null` if the description names no segment); `null` without references |
| `settings` | object | the options, see below |
| `length_bin` | integer | width of the bins of the length histograms (10 bases) |
| `length_max` | integer | the last bin of a histogram contains all lengths from this value on (4000) |
| `short_fraction` | number | a complete amplicon is short if shorter than this fraction of its reference (0.9) |
| `all` | section | statistics of all reads |
| `groups` | array | for each group of input files (a file given directly, a directory, or each barcode directory of a `fastq_pass` directory): `name`, `files` (number of input files), `processed` (`false`: not processed yet, then no further keys), and the keys of a section |

Settings:

| key | type | description |
|---|---|---|
| `max_errors` | object | the errors allowed in a primer (`-e`): `mode` (`"fixed"`), `value` (below 1: a fraction of the length of the primer, otherwise a number of bases) |
| `window` | integer | `-w`: a primer must lie within this many bases of its end of the read |
| `min_length` | integer | `-m`: shorter reads after the trimming are not written |
| `probe_length`, `probe_errors`, `terminal_window`, `exhaustive_search` | integer, number, integer, boolean | the options of the classification against the references (used only with `-r`) |
| `tags` | boolean | the tags are appended to the titles of the reads (`false`: `--no-tags`) |
| `threads` | integer | `-j` |
| `compress_level` | integer | `--compress-level`: the gzip compression level of the outputs |

Section (statistics of a set of reads):

| key | type | description |
|---|---|---|
| `reads` | object | `in`, `bases_in` (reads read and their bases), `written`, `bases_written`, `both_primers`, `one_primer`, `no_primer` (reads by number of ends at which a primer was found), `dimers` (both primers and nothing between them), `too_short` (not written: empty or shorter than `-m` after the trimming; dimers included), `minus_strand` (reads whose primers were found on their reverse complement), `with_internal_match` (reads with primer matches inside them), `adapter_start`, `adapter_end` (reads with a Nanopore adapter removed at their start, at their end: searched only at the ends without a primer, `--adapters`) |
| `primers` | object | for each end (`left`, `right`): the reads in which each primer of the scheme was found at that end, by primer name |
| `ends` | object or null | with references: for each end, the ends found by origin: `terminal`, `internal`, `unknown` (not located in the references); `null` without references |
| `internal_sites` | array or null | with references: the internal priming sites, by decreasing number of reads (at most 1000): `reference`, `end` (`left` or `right`), `position` (1-based position on the plus strand of the reference of the first base after the primer, as in the `ml`/`mr` tags of the reads, where most read ends of the site mapped), `from`, `to` (the range of the positions of the site: positions at most 5 bases apart are one site), `reads`; `null` without references |
| `short_complete` | integer or null | with references: the complete amplicons shorter than 90% of their reference (as in `summary.tsv`); `null` without references |
| `complete_amplicons` | array | the complete amplicons (reads with the primers of both ends) of each group of the lengths with at least one amplicon: `group` (the reference to which the bases after the left primer, or else the right primer, mapped; `not_located`; or `all` without references), `reference_length` (`null` if not a reference), `amplicons`, `main_length` (start of the most frequent bin), `short` (amplicons shorter than 90% of the reference; `null` if not a reference), `short_main_length` (start of the most frequent bin of the short amplicons; `null` if none), `histogram` (the counts of the bins of the length between the primers, from 0, up to the last non-empty bin) |

For example, the internal priming sites of each group with `jq`:

    jq -r '.groups[] | .name as $g | .internal_sites[]? |
           [$g, .reference, .end, .position, .reads] | @tsv' \
        out/flucut_summary.json
