# JSON summary

flucut always writes a machine-readable summary of the run to
`outdir/flucut_summary.json`, with the same information as the HTML
report (`outdir/flucut_report.html`) and the dashboard
(`outdir/flucut_dashboard.html`). The file is written to a temporary
file and then renamed, so that it can be read at any time.

## Format (version 1)

Top level:

| key | type | description |
|---|---|---|
| `format` | integer | version of this format (1) |
| `program`, `version` | string | `"flucut"` and its version |
| `updated` | string | time of writing (ISO 8601) |
| `running` | boolean | `true` while the run is in progress (the summary is rewritten after each group of input files, in the live mode after each file); `false` in the summary of the end of the run |
| `command` | string | the command line |
| `scheme` | object | the primer scheme: `name`, and the names of the primers of the `left` and `right` ends of an amplicon, arrays of strings, in the order of the scheme |
| `trim_policy` | string | `"primer"` (no references) or `"terminus"` (with `-r`: conserved termini kept at terminal ends) |
| `references` | object or null | with `-r`: `file` (the path given), `sequences` (array of `name`, `length`, `virus` (`"A"` or `"B"`, from the FASTA description), `segment` (1-8, from the description and read in the numbering of that type, see the manual) and `gene` (e.g. `HA`), `null` if the description names no segment); `null` without references |
| `settings` | object | the options, see below |
| `length_bin` | integer | width of the bins of the length histograms (10 bases) |
| `length_max` | integer | the last bin of a histogram contains all lengths from this value on (4000) |
| `short_fraction` | number | a complete amplicon is short if shorter than this fraction of its reference (0.9) |
| `min_depth` | integer | the depth from which a position is well covered (50) |
| `complete_fraction` | number | a segment is complete if at least this fraction of its positions is well covered (0.9) |
| `second_subtype_min_excess`, `second_subtype_min_share` | number | the reads of a second HA or NA subtype are above the threshold if they exceed those expected by chance by at least this number of reads (2.5: 3 reads where none is expected) and by at least this fraction (0.01) of the reads counted (see `subtype_reads` below) |
| `all` | section | statistics of all reads |
| `groups` | array | for each group of input files (a file given directly, a directory, or each barcode directory of a `fastq_pass` directory): `name`, `files` (number of input files), `processed` (`false`: not processed yet, then no further keys), and the keys of a section |

Settings:

| key | type | description |
|---|---|---|
| `max_errors` | object | the errors allowed in a primer (`-e`): `mode` (`"fixed"`), `value` (below 1: a fraction of the length of the primer, otherwise a number of bases) |
| `platform` | string | `-p`: `"nanopore"` or `"illumina"` |
| `window` | integer | `-w`: a primer must lie within this many bases of its end of the read |
| `min_length` | integer | `-m`: shorter reads after the trimming are not written |
| `probe_length`, `probe_errors`, `terminal_window`, `exhaustive_search` | integer, number, integer, boolean | the options of the classification against the references (used only with `-r`) |
| `segment` | integer or null | `--segment`: the segment (1-8) whose references alone were used for the classification; `null`: all segments |
| `report_segment` | integer or null | `--report-segment` (or `--segment`): the segment (1-8) shown in the per-reference parts of the HTML report; `null`: all. The JSON summary itself always describes all the references |
| `tags` | boolean | the tags are appended to the titles of the reads (`false`: `--no-tags`) |
| `split` | boolean | the reads which contain several amplicons are split (`false`: `--no-split`) |
| `cut_deletions` | boolean | `--cut-deletions`: a read with an internal deletion is written as two records, its two flanks, cut at the junction |
| `by_category` | boolean | `--by-category`: the reads of each input are written to five files, one for each category of reads (full templates, internal priming, internal deletions, complete amplicons placed on no segment, other) |
| `threads` | integer | `-j` |
| `compress_level` | integer | `--compress-level`: the gzip compression level of the outputs |

Section (statistics of a set of reads):

| key | type | description |
|---|---|---|
| `reads` | object | `in`, `bases_in` (reads read and their bases), `written`, `bases_written`, `both_primers`, `one_primer`, `no_primer` (reads by number of ends at which a primer was found), `dimers` (both primers and nothing between them), `too_short` (not written: empty or shorter than `-m` after the trimming; dimers included), `minus_strand` (reads whose primers were found on their reverse complement), `with_internal_match` (reads with primer matches inside them), `split` (reads split into several amplicons: `in` and `bases_in` count the reads of the input, the other counts each part as a read), `cut_at_junction` (reads with an internal deletion cut at their junction, `--cut-deletions`: `written`, `bases_written`, `too_short` and `written_by_category` count each flank as a read), `adapter_start`, `adapter_end` (reads with a Nanopore adapter removed at their start, at their end: searched only at the ends without a primer, `--adapters`) |
| `read_lengths` | object | for the reads of the `input` and for the reads `written`: `reads`, `bases`, `min`, `max`, `mean`, `median` and `n50` (the length such that half of the bases are in reads at least that long); the median and the N50 are exact to 20 bases |
| `read_quality` | object | for `input` and `written`: `reads` (those with qualities), `median` (of the qualities of the reads, to 1 unit) and `q10`, `q20` (the fractions of the reads with a quality of at least 10, 20); the quality of a read is the Phred score of the mean error probability of its bases; `null` for reads without qualities |
| `sequencing` | object | what the titles of the reads say about the run: `flow_cell`, `run`, `basecall_model`, `sample`, `experiment` (the fields of MinKNOW), `instrument` (Illumina) and `start_time` (the earliest among the titles read); each a string, `"several"` if the reads differ, or `null` if not known. The title of one read in a thousand is read |
| `sample` | object or null | with references, the sample itself: `state` (`complete`: all segments expected are complete; `partial`; `incomplete`: some are missing; `few_reads`, `not_located`: not complete, and it cannot be judged, with fewer than 1000 reads or half or more of the read ends not located), `virus` (`"A"` or `"B"`: the type of the references its reads were located in), `subtype` (string or null: the subtype for influenza A, e.g. `"H5N1"`, and `"B"` with its lineage where a reference names one, e.g. `"B/Victoria"`, for influenza B), `segments_expected`, `segments_found`, `segments_complete`, `ends_not_located` (fraction), `reads_placed` (fraction of the reads, after splitting, placed on a segment of the references), `subtype_reads` (the reads of each HA and NA subtype, see below), `flags` (what deserves attention: `few_reads`, `no_primer`, `not_located`, `few_influenza_reads`, `internal_priming`, `internal_deletions`, `cleavage_site`, `second_ha_subtype`, `second_na_subtype`) and `segments`: for each segment expected, `segment` (its number in the type of the sample: influenza B numbers PB1 1 and PB2 2), `gene`, `state` (`missing`, `low`, `partial`, `complete`) and, if found, `reference`, `length`, `mean_depth`, `median_depth`, `median_full_depth` (the median depth counting only the full templates, the reads whose amplicon started at a segment end), `covered` and `well_covered` (fractions of the positions with a read, with at least `min_depth` reads), `internal`, `deletion` and `minus_strand` (fractions of the depth), `complete_amplicons`, `short_complete`. Depths are estimated from the positions of the reads, not from a mapping. `null` without references |
| `cleavage_site` | object or null | with references and an HA reference containing the fusion peptide: the cleavage site of the haemagglutinin as the reads carry it (see `docs/cleavage-site.md`): `chemistry` (`"Illumina"`, `"R10"`, `"R9"` or `"unknown"`), `reads_searched`, `reads` (the reads the site was read from), `distinct_motifs` and `distinct_sites` (the distinct spellings and the distinct residues among those reads, the readings of single misread reads included, so large in nearly every sample), `supported_sites` (the residues carried by at least 10% of the reads and at least 2 reads: the sites the sample holds), `motif` (the 24 bases, the most frequent spelling) and `protein` (the 8 residues), `motif_reads`, `support` (fraction of `reads` carrying those residues), `multibasic_share` (fraction of the reads whose own site is multibasic), `basic_residues`, `control_reads` and `control_support` (the same reads in two windows beside the site, where the sequence is ordinary), `reads_disagree_beside_the_site`, `fusion_peptide` with `fusion_peptide_reads` and `fusion_peptide_support` (the anchor as the reads carry it) and `reference_fusion_peptide`, `reads_as` (`multibasic`, `monobasic` or `unreadable`), `reference_site` (the site of the HA reference of the sample) and `disagrees_with_reference`, `verdict` (`firm`, `confirm` or `too few reads`), and `sites`: the most frequent sites in the same reads, each `protein`, `motif`, `reads`, `share` and `reads_as`. `null` where the site was not read |
| `primers` | object | for each end (`left`, `right`): the reads in which each primer of the scheme was found at that end, by primer name |
| `written_by_category` | object or null | with references: the reads written, by category (the `ct` tag): `full` (full templates), `internal` (internal priming), `deletion` (internal deletions), `unplaced` (a complete amplicon the references place on no segment), `other` (no 5' primer, or not located); their sum is `reads.written`; `null` without references |
| `ends` | object or null | with references: for each end, the ends found by origin: `terminal`, `internal`, `unknown` (not located in the references); `null` without references |
| `internal_sites` | array or null | with references: the internal priming sites, by decreasing number of reads (at most 1000): `reference`, `end` (`left` or `right`), `position` (1-based position on the plus strand of the reference of the first base after the primer, as in the `ml`/`mr` tags of the reads, where most read ends of the site mapped), `from`, `to` (the range of the positions of the site: positions at most 5 bases apart are one site), `reads`; `null` without references |
| `short_complete` | integer or null | with references: the complete amplicons shorter than 90% of their reference (as in `summary.tsv`); `null` without references |
| `complete_amplicons` | array | the complete amplicons (reads with the primers of both ends) of each group of the lengths with at least one amplicon: `group` (the reference to which the bases after the left primer, or else the right primer, mapped; `not_located`; or `all` without references), `reference_length` (`null` if not a reference), `amplicons`, `main_length` (start of the most frequent bin), `short` (amplicons shorter than 90% of the reference; `null` if not a reference), `short_main_length` (start of the most frequent bin of the short amplicons; `null` if none), `histogram` (the counts of the bins of the length between the primers, from 0, up to the last non-empty bin) |
| `coverage` | array or null | with references: for each reference with located reads, the coverage by read origin, as drawn in the report: `reference`, `length`, `positions` (up to 200 positions, 0-based, evenly spaced along the reference) and `depth`, with the depth at those positions of each class of reads with reads on the reference: `terminal_plus`, `terminal_minus` (full templates, by strand), `internal_plus`, `internal_minus` (internally primed), `no_primer`, `not_located` (the two kinds of other reads) and `deletion` (internal deletions). It is estimated from the positions of the reads, not from a mapping; `null` without references |

Reads of each HA and NA subtype (`sample.subtype_reads`): an object with
the keys `ha` and `na`, each an object or `null`; `null` as a whole
where nothing was counted (influenza B, references without subtypes or
with one reference per segment, and `all` in a run of several samples,
whose subtypes may differ). The counts are those of the reads from
which the references of the sample were chosen: the first 2000 reads of
its first file (up to 20,000 if a segment has too few probes located
in them).

| key | type | description |
|---|---|---|
| `leading` | string | the subtype with most probes located in its references, e.g. `"H5"` (`"N1"` for the NA) |
| `second` | string or null | the other subtype whose reads exceed most those expected by chance; `null` if no read is of another subtype |
| `reads` | integer | the reads with several voting probes in the gene (a probe votes if it is located in the references and no reference of another subtype matches it as well) |
| `leading_reads`, `second_reads` | integer | of them, the reads whose probes all vote for the leading subtype, for the second one |
| `second_expected` | number | the reads expected to be all for the second subtype if each probe voted independently of its read (a read with *k* voting probes is all for a subtype with the fraction *p* of the votes with probability *p*^*k*) |
| `second_excess` | number | `second_reads` minus `second_expected` |
| `second_share` | number | `second_excess` as a fraction of `reads` |
| `above_threshold` | boolean | `second_excess` is at least `second_subtype_min_excess` and `second_share` at least `second_subtype_min_share`; the sample then has the flag `second_ha_subtype` or `second_na_subtype` |

Reads of a second subtype can come from a mixed infection, from a
contamination of the sample, or from barcode crosstalk between the
samples of a run; the summary says what was observed, not which of them
it is. Two viruses of one subtype are not told apart, and the other six
segments are not looked at.

For example, the internal priming sites of each group with `jq`:

    jq -r '.groups[] | .name as $g | .internal_sites[]? |
           [$g, .reference, .end, .position, .reads] | @tsv' \
        out/flucut_summary.json

The state, the subtype and the segments at a depth of 50 or more of
each sample:

    jq -r '.groups[] | [.name, .sample.state, .sample.subtype,
           (.sample.segments | map(.well_covered // 0 | . * 100 | floor)
            | join(" "))] | @tsv' out/flucut_summary.json

The samples with reads of a second HA subtype above the threshold:

    jq -r '.groups[] | select(.sample.subtype_reads.ha.above_threshold)
           | [.name, .sample.subtype_reads.ha.leading,
              .sample.subtype_reads.ha.second,
              .sample.subtype_reads.ha.second_reads,
              .sample.subtype_reads.ha.reads] | @tsv' \
        out/flucut_summary.json
