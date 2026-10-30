# flucut real-data search — Nanopore only

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

## Step 4: Nanopore Uni12/Uni13 reads checked directly — artefact FOUND in reads

Protocol screen (SRA experiment XML, ENA) of the H5 Nanopore amplicon
projects: **PRJNA1406649** (Cornell University, "Oxford Nanopore
sequencing of Highly Pathogenic Avian Influenza H5N1 reporter viruses",
clade 2.3.4.4b; design: "MBTUni Primers were used to Amplify all 8
segments", SQK-LSK109; 6 runs, SRR39646234-SRR39646239) uses the
Uni12/Uni13 (MBTuni) multisegment RT-PCR. Others looked at: PRJNA1312032
(Penn State, retail milk, no primer stated), PRJEB79943 (Bern, mRT-PCR,
2 H5N1 runs, not yet checked), PRJNA1173216 (USDA, SISPA, not Uni12).

Two smallest Cornell runs were downloaded (ENA) and the HA site counted
in the raw reads (scripts/count_codon_nanopore.py, exact match on the 12
template bases downstream of the Uni12-like site, both strands; the
site is HA nt 1050 of PX279898.1, codon nt 1054-1056). Reads are split
by whether the site lies at the read start (<= 60 bases: internal
priming) or inside the read (the read spans it: the template).

| run | site at read start: AAA / ATA | read spans site: AAA / ATA |
|---|---|---|
| SRR39646237 (12,532 reads) | 512 / 5 | 9 / 837 |
| SRR39646239 (15,482 reads) | 122 / 4 | 3 / 508 |

Reads spanning the site carry the template codon ATA (Ile); reads that
start at the site carry AAA (Lys), the primer's own bases. This is the
artefact in raw Nanopore reads. flucut (`-r isolate.fasta`, default
`--trim primer`) reports an HA length peak at ~600 bases (internally
primed amplicons from the site to the 3' end) and 1661 / 490
left-internal reads (summary.tsv).

Caveat: the reference used for flucut and the count is PX279898.1
(H5N1 cattle isolate, carries AAA); the isolates of the Cornell runs
are reporter viruses (NS reporter) of the same clade, so the site
flank matches exactly. The effect on a consensus (IRMA untrimmed vs
trimmed) was not run here (no IRMA in this container).

**Candidate benchmark data: SRR39646237 and SRR39646239 (PRJNA1406649,
Nanopore, H5N1, Uni12/Uni13).**

## Step 5: natural (non-engineered) isolates in Nanopore reads

Searched for natural H5N1 isolates with the artefact in reads, same
counter (scripts/count_codon_nanopore.py, reference PX279898.1, site
HA nt 1050).

**Positive: PRJNA1114404** (Cornell University, 41 runs, GridION,
"targeted whole genome sequencing of Influenza A virus H5N1", "PCR
amplification using MBTuni primers", host samples 2024, spillover of
HPAI H5N1 to dairy cattle). The five smallest runs, all positive:

| run | BioSample | sample | site at read start: AAA / ATA | read spans site: AAA / ATA |
|---|---|---|---|---|
| SRR29165103 | SAMN41527190 | 078707-24-14 | 456 / 15 | 5 / 956 |
| SRR29165110 | SAMN41527219 | 067041-24-1 | 175 / 2 | 0 / 153 |
| SRR29165113 | SAMN41527216 | 063224-24-10 | 264 / 7 | 0 / 552 |
| SRR29165115 | SAMN41527215 | 063224-24-1 | 372 / 4 | 3 / 824 |
| SRR29165125 | SAMN41527188 | 078707-24-12 | 120 / 3 | 2 / 320 |

(out/cornell2024_codon_counts.tsv.) None of these BioSamples has a
motif-covering consensus among the deposited records screened in step
1-2, so the artefact is in the reads and no K-consensus is public for
them. The other 35 runs of the project were not scanned.

**Negative:** PRJNA1312032 (Penn State, retail milk H5N1, 3 runs): reads
spanning the site are all ATA (SRR35186256: 69; SRR35186257, first
150 MB: 1153), no read starts at the site; SRR35186258 (first 150 MB)
has no read with the site flank. No Uni12 internal-priming signal.
Not yet scanned: PRJNA1400594 (Michigan, 63 wild-bird H5N1 runs,
"universal IAV primers", avgLength 427), PRJEB79943 (Bern, 2 H5N1 runs).

## Files
- out/nanopore_runs.tsv — merged SRA runinfo (4123 runs)
- out/motif_hits.tsv — all 243 motif-covering consensus records + residue
- scripts/ — all scripts (runinfo, elink, esearch, efetch, screen)
- out/cornell_codon_counts.tsv — codon counts of step 4
- out/cornell_PRJNA1406649_ena.tsv — ENA run list of PRJNA1406649
- out/cornell2024_codon_counts.tsv — codon counts of step 5
- dl/ — downloaded data (gitignored)
