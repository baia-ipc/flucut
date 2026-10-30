# flucut real-data search — Nanopore only (in progress)

All queries, downloads and runs recorded. GISAID not used (not accessible).
Date of work: 2026-10-30. Attribution: H5N1 HA "I348K" (GLFGAIAG→GLFGAKAG)
Uni12 internal-priming artefact; see docs/internal-priming.md.

## Method

The artefact is a K at the 6th residue of the HA fusion-peptide motif
`GLFGA[I]AG` (highly conserved across influenza A HA2). The variable codon
`ATA` (Ile) is overwritten by the Uni12 primer's `AAA` (Lys) in internally
primed amplicons. Screen = 6-frame translation, regex `GLFGA(.)AG`,
report the 6th residue. Positive control **PX279898.1** (known Illumina
artefact consensus, A/cattle/CA/25-020784-009-original/2025) correctly
returns `GLFGAKAG` at segment nt 1054 — screen validated.

## Step 1–2: Nanopore SRA runs → their deposited consensus → screen

Query (SRA, 2026-10-30):
`"Influenza A virus"[Organism] AND "OXFORD_NANOPORE"[Platform]`
→ 4123 runs (all Nanopore), 3240 BioSamples, 65 BioProjects.
Library: WGS 1553, FL-cDNA 1092, OTHER 806, AMPLICON 621; selection
RT-PCR 3096, other 807, PCR 130. (runinfo in out/nanopore_runs.tsv.)

Linked the 3240 BioSamples to GenBank (elink biosample→nuccore):
1044 BioSamples carried deposited sequences; 7433 nuccore UIDs,
7412 FASTA fetched. Screened all 7412.

**Result: 243 records cover the fusion-peptide motif; ALL 243 are Ile
(`GLFGAIAG`). 0 carry the `GLFGAKAG` artefact.**

So **no GenBank consensus deposited from a public Nanopore influenza A run
carries the deposited I348K artefact.** This is consistent with the earlier
Illumina scan (benchmark/RESULTS.md): the 8 artefact-linked runs there were
all Illumina (USDA dairy-cattle surveillance, PRJNA1102327). The deposited
K-consensus artefact is, in the public data, Illumina-associated.

(Caveat: this screens *deposited consensus*. A Nanopore run can still carry
the artefact *in its reads* while its deposited/QC'd consensus reads Ile,
or while no consensus was deposited. That is Step 3–4, below, and is the
point of checking reads directly with flucut.)

## Step 4 (next): Nanopore Uni12/Uni13 multisegment RT-PCR read datasets

3232 of the 4123 Nanopore runs are AMPLICON or RT-PCR/PCR selection —
candidate multisegment RT-PCR data to check directly with flucut even
without a deposited K consensus. Largest projects (out/nanopore_runs.tsv):

| runs | BioProject | center |
|---|---|---|
| 1068 | PRJEB56915 | University of Oxford |
| 998 | PRJNA1213174 | Hong Kong Polytechnic University |
| 200 | PRJNA1456957 | Public Health Agency of Canada |
| 144 | PRJEB79943 | Inst. for Infectious Diseases, Bern |
| 61 | PRJNA1153990 | CDC-FLU |
| 41 | PRJNA1114404 | Cornell University |

To do: identify which use Uni12/Uni13 (MBTuni) multisegment RT-PCR, pick
1–2 small H5N1 runs, run `flucut -r <refs>`, and count reads starting at
the site (primer `AAA`) vs reads spanning it (`ATA`).

## Files
- out/nanopore_runs.tsv — merged SRA runinfo (4123 runs)
- out/motif_hits.tsv — all 243 motif-covering consensus records + residue
- scripts/ — all scripts (runinfo, elink, esearch, efetch, screen)
- dl/ — downloaded data (gitignored)
