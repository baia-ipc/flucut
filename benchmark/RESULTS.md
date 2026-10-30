# Results

## 2026-10-30: IRMA checks (simulated reads)

What was run: `benchmark/run-irma-checks.sh` (see `README.md`), on the
genome of `make-genome.sh` (downloaded 2026-10-30), with the defaults of
`simulate.py` (seed 1, 400 reads per segment, 3200 reads, 239 of them
started at the internal site of HA).

Versions: flucut at commit 4eec29e (with the benchmark scripts of this
commit; the program is unchanged); IRMA 1.4.0 (commit 1311c81),
IRMA-core 0.11.0, IRMA-viz 0.2.2, LABEL 0.7.2; module FLU-minion.

Runs: `raw` (untrimmed), `flucut` (default trimming), `notags` (the same
with `--no-tags`), `terminus` (`--trim terminus -r` the genome); each
with FLU-minion as it is and with `SKIP_E=0` (`_elong`).

### The I348K site

The codon of the site in the plain (majority) consensus of HA:

| run | codon |
|---|---|
| raw, raw_elong | **AAA** (Lys: the artefact) |
| flucut, notags, terminus (both modes) | ATA (Ile: the truth) |

Untrimmed, at the middle base of the codon 61% of the reads carry the A
of the primer and 37% the T of the template (IRMA allele table of the
raw run); the majority consensus calls AAA, the amended consensus an
ambiguity (W). With flucut the site is the truth in all runs. With 0.6
of the HA reads internally primed (a simulation parameter, not a
measured share), the untrimmed reads reproduce the published artefact.

### Tags in the read titles

`flucut` and `notags` give the same IRMA read counts and the same
amended consensus in all 8 segments, with both configurations: IRMA
accepts the tags. `--no-tags` stays an option, not the default.

### The segment ends

Bases of the true segment missing (or beyond it) at each end of the
amended consensus, and errors:

| run | 5' end | 3' end | errors |
|---|---|---|---|
| any input, FLU-minion as it is | 20-45 missing | 20-58 missing | none (raw: 2 ambiguous at the HA site) |
| raw_elong | 0-28 **beyond** the truth | 13-32 beyond | 1-6 insertions per segment, ambiguous bases |
| flucut_elong | 12 missing | 13 missing | none (segment 3: one N appended at the 3' end) |
| terminus_elong | complete | complete | 1 ambiguous base per segment; segment 2: 1 substitution; segments 3 and 4: one N appended |

- With FLU-minion as it is, IRMA does not extend the consensus beyond
  its references: the same 20-58 bases are missing at the ends whatever
  the trimming, so keeping the termini changes nothing.
- With the elongation, the consensus of the flucut reads ends exactly
  where the primers began (12 and 13 bases short), without errors; the
  untrimmed reads extend it into flanks and primer tails, with
  insertions.
- Keeping the termini (`--trim terminus`) gives complete ends, but the
  bases there are the primer's: the 4th base is an ambiguity (R) in all
  segments, from the degenerate Uni12 primer (the even mix of A and G
  is a simulation assumption), and in segment 2 the 3' terminus of PR8,
  `CCTTGTTCCTACT`, is overwritten by the primer, `CCTTGTTTCTACT` (real
  sequences, no assumption). This supports removing every primer by
  default.

### Not yet measured

The other tools of `TOOLS.md`, the read-level measures (primer bases
left, template bases removed, reads lost: the truth is in
`truth.tsv.gz`), other shares of internal priming, other error rates,
and real data.
