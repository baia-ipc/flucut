# Real-data benchmark: Uni12 internal-priming artefact in public H5N1 data

**Date:** 2026-10-30  
**Tool:** flucut (commit 3bfb35c)

## Goal

Find public influenza sequencing data where the Uni12 internal-priming
artefact produces a I→K substitution in the HA fusion peptide
(GLFGAIAG → GLFGAKAG), link to the raw reads in SRA, and verify with
flucut.

## Coordinates (reference: KF918483.1, H5N1 HA)

| Numbering        | Position of I/K residue |
|-------------------|------------------------|
| Full ORF (Met=1)  | aa 352                 |
| Mature HA (sp=16) | aa 336                 |
| HA2 (after RRKKR) | aa 6                   |
| Segment nt        | 1054–1056 (ATA=Ile wt) |

The artefact codon is **AAA** (Lys) instead of **ATA** (Ile).

## Scan of 29,997 public H5N1 HA consensus sequences

Source: NCBI nuccore, query `"Influenza A virus"[Organism] AND H5N1
AND (segment 4 OR hemagglutinin) AND 1600:2000[SLEN]`, downloaded
2026-10-30.

| Category       | Count  |
|----------------|--------|
| Wildtype GLFGAIAG | 29,777 |
| **GLFGAKAG (K artefact)** | **109** |
| Other GLFGA.AG variants | 27 |
| No motif found | 84 |

The 109 GLFGAKAG records span multiple countries and years (2006–2026),
including cattle (US), poultry (Bangladesh, Myanmar), wild birds
(US, Canada, Germany, South Georgia, Antarctica, Brazil, Namibia,
Yemen, Switzerland), mammals (mink/Finland, fur seal/Antarctica), and
one human case (Nanjing 2015). Full list in `mut_accessions.txt`.

The 27 "other" variants include GLFGAXAG (ambiguous base calls) and
a few GLFGAVAG records; these likely also reflect sequencing or
assembly issues in the fusion peptide region.

## SRA linkage

Of the 109 GLFGAKAG records, 8 link to public SRA runs via
BioSample, all from the USDA dairy-cattle H5N1 surveillance project
(PRJNA1102327):

| SRA Run       | Isolate                                    | Spots   |
|---------------|--------------------------------------------|---------|
| SRR31168528   | A/cattle/CA/24-028986-002-original/2024    | 460,456 |
| SRR34913296   | A/cattle/CA/25-020784-009-original/2025    | 282,828 |
| SRR37221769   | A/cattle/CA/26G01054-008-original/2024     | 413,846 |
| SRR37221768   | A/cattle/CA/26G01054-010-original/2024     | 577,297 |
| SRR39413549   | A/cattle/ID/26G08404-029-original/2026     | 616,830 |
| SRR39565236   | A/cattle/UT/26G08463-001-rep2/2026         | 371,031 |
| SRR40141856   | A/cattle/UT/26G10577-001-original/2026     | 492,680 |
| SRR40141855   | A/cattle/UT/26G10578-001-original/2026     | 660,250 |

The remaining 101 records have no linked SRA run (consensus-only
submissions).

## flucut analysis of SRR34913296

**Run:** SRR34913296 (A/cattle/CA/25-020784-009-original/2025, paired
Illumina NextSeq 2000, 282,828 spots).

**References:** The isolate's own 8-segment consensus (PX279895–
PX279902).

**Command:**
```
flucut -o flucut_out -r isolate_all8.fasta R1.fastq.gz R2.fastq.gz
```

### Summary (both read files combined)

| Metric                          | R1       | R2       |
|---------------------------------|----------|----------|
| Reads in                        | 282,828  | 282,828  |
| Reads out                       | 282,821  | 282,821  |
| Both primers found              | 186      | 216      |
| One primer found                | 43,345   | 41,529   |
| No primer found                 | 239,297  | 241,083  |
| Left ends: internal priming     | 5,819    | 5,961    |
| Right ends: internal priming    | 269      | 247      |

### HA segment internal-priming hotspots (PX279898.1)

1,779 reads mapped to internal-priming sites on HA. Top positions
(segment coordinates, 1-based):

| Position | Reads | Offset from K codon (nt 1054) |
|----------|-------|-------------------------------|
| 1130     | 1,312 | +76 nt                        |
| 1062     | 272   | +8 nt                         |
| 1539     | 76    | +485 nt                       |
| 1267     | 46    | +213 nt                       |
| 1085     | 42    | +31 nt                        |

The **position 1062** hotspot is immediately downstream of the AAA
(Lys) codon at nt 1054–1056 — consistent with the Uni12 primer
binding internally at the conserved AGCAAAAGCAGG-like motif in the
fusion peptide coding region, producing truncated amplicons that
start near nt 1062 instead of at the segment terminus.

Position 1130 (the largest hotspot) is 76 nt further downstream,
representing a second internal-priming site in the same region.

### Internal priming across all 8 segments

flucut detected internal priming on every segment, with the highest
counts on:

- **Segment 7 (M):** 1,712 reads at position 892
- **Segment 4 (HA):** 1,312 reads at position 1130, 272 at 1062
- **Segment 5 (NP):** 849 reads at position 828
- **Segment 8 (NS):** 842 reads at position 489

## Conclusions

1. **109 of ~30,000 public H5N1 HA consensus sequences carry the
   GLFGAKAG motif** (0.36%), spanning multiple hosts, countries, and
   years.

2. **8 of these 109 link to public SRA runs** (all USDA dairy cattle
   surveillance). The other 101 are consensus-only deposits with no
   raw reads in SRA.

3. **flucut confirms internal priming on HA near the fusion peptide
   region** in SRR34913296, with a hotspot at nt 1062 — 8 bases
   downstream of the I→K codon. This is direct evidence that Uni12
   can prime internally at the GLFGA motif, producing truncated
   reads that bias the consensus toward K (AAA) over I (ATA).

4. The artefact is not limited to HA: flucut finds internal-priming
   hotspots on all 8 segments, consistent with the Uni12 primer
   matching conserved motifs throughout the genome.

## Files

| File | Description |
|------|-------------|
| `KF918483.fasta`, `KF918483.gb` | H5N1 HA reference |
| `cambodia_h5n1_ha.fasta` | 137 Cambodia H5N1 HA sequences (all wt) |
| `h5n1_ha_all.fasta` | 29,997 H5N1 HA sequences |
| `scan_results.txt` | Translation scan output |
| `mut_accessions.txt` | 109 GLFGAKAG accessions |
| `other_accessions.txt` | 27 other-variant accessions |
| `sra_runs.tsv` | 8 SRA runs linked to GLFGAKAG records |
| `refs_all8.fasta` | RefSeq 8-segment references |
| `isolate_SRR34913296_all8.fasta` | Isolate-specific 8 segments |
| `sra_download/SRR34913296_{1,2}.fastq.gz` | Raw reads |
| `flucut_out_SRR34913296_v2/` | flucut output with isolate refs |