# The HA cleavage site, read from the reads

The haemagglutinin is made as one protein and cut in two, HA1 and HA2,
before the virus can fuse with a cell. The residues at the cut are the
cleavage site. In the highly pathogenic avian H5 and H7 viruses it is a
run of basic residues, for example `PQRERRRKR`; in the low pathogenic
ones a single one, for example `PQRETR`. That difference is what
distinguishes the two pathotypes, so a consensus which is wrong there
reports the wrong virus.

The run of basic residues is encoded by an A/G-rich stretch, which is
exactly the kind of low-complexity sequence that sequencing reads get
wrong, by an insertion or a deletion. An assembler called from few reads
carries such an error into the consensus; the reading frames of HA1 and
HA2 then differ and the site is nonsense, or, with three bases too many,
the run of basic residues is silently one residue longer.

flucut reports what the **reads** carry at the site, with how many of
them agree. This is not an assembly and does not replace one: it says
what the data supports before any consensus is built, and whether the
consensus that follows can be trusted there.

## Where it is read from

HA2 begins directly after the cleavage site with the fusion peptide,
which every haemagglutinin has and which lies outside the low-complexity
stretch. It is `GLFGAIAG` in most of them, `GIFGAIAG` in the H3 of
horses and pigs, `GLLGAIAG` in some others and `GFFGAIAG` in influenza
B, so flucut looks for the pattern `G..GAIAG`, in which the two
residues that vary are free. The pattern occurs exactly once in each of
the 89 HA references of the reference sets flucut ships, of every
subtype H1 to H16 and of influenza B. flucut takes those 24 bases of
the HA reference as an **anchor**, finds the anchor inside each read
(edlib, infix, both orientations, the same error allowance as the
probes), and reads the 24 bases immediately before the match, as the
read carries them.

Reading beside an anchor, rather than from a mapping, has two
consequences worth knowing:

- the frame is set by the anchor, which abuts the site, so an insertion
  or a deletion anywhere else in the read cannot shift what is read;
- nothing is aligned to anything: one search of a 24-base pattern per
  read, no assembly, so the result is available while the run is still
  sequencing.

The reads are then gathered by the **residues** they carry, not by the
bases. The residues are what decides the pathotype, and the reads agree
on them far more often than on the spelling, which differs by a base
wherever a read has an error.

## Which reads are read

What reaches the site differs between the platforms, so the reads taken
differ too.

| platform | reads read | why |
|---|---|---|
| Nanopore (`-p nanopore`) | the full templates located on the HA reference | a read of a whole amplicon reaches the site only when its amplicon started at a segment end: the internally primed reads begin just after the site (the Uni12-like site lies a few bases after it) and the genomes with an internal deletion have removed it |
| Illumina (`-p illumina`) | every read | a read is a fragment of an amplicon, so a read which covers the site carries no primer, has no category and is not located at all; the anchor itself says which reads cover the site |

## How far it can be trusted

Two things are measured from the same reads, and both have to hold.

**Do the reads agree at the site?** The share of them carrying the
majority residues is the **support**.

**Do they agree where nothing is difficult?** Beside the anchor, on the
HA2 side, lie two windows of the same width as the site, taken from the
same reads. HA2 there is ordinary coding sequence. The share of reads
agreeing in them is the **control support**, and it says what this
sample's reads agree on when the sequence is not against them.

The second is what tells a site the reads cannot resolve from reads
that agree nowhere. A sample whose control is high and whose site is low
has a genuinely hard site. A sample low in both has reads too noisy, or
holding more than one virus, for a claim about any position.

| | needed to be firm |
|---|---|
| reads over the site | 20 |
| support at the site | 0.45 |
| support in the control windows | 0.60 |

The verdict is one of

- **firm**: enough reads, enough of them agree at the site, they agree
  beside it, and they read like the reference of the sample;
- **confirm**: one of those fails. The residues are still reported, with
  the reason. Confirm before reporting a pathotype.
- **too few reads**: fewer than 20 reads reach the site; nothing is said
  about it.

Neither threshold names a platform or a chemistry. The chemistry is
still read from the basecalling model in the read titles and reported,
because it is worth knowing, but it decides nothing.

Two further numbers are reported next to the majority, and are often
more useful than it:

- the share of the reads whose **own** site is multibasic, whatever its
  spelling. The reads can disagree on the exact residues and still agree
  that the site is multibasic, which is what the pathotype turns on.
- the **other sites** found in the same reads. The fusion peptide is the
  same in every subtype, so the anchor finds the site of every
  haemagglutinin present in the sample. A barcode which holds the reads
  of more than one of them shows here as two sites with substantial
  shares, and when the majority does not read like the HA reference
  chosen for that sample, flucut says so and never calls the result
  firm.

## Where the numbers come from

Measured on 104 Nanopore samples of avian H5N1, H5N2 and H9N2 (94 R10,
10 R9) of one laboratory and one protocol, run through flucut itself.
Twelve of them have an IRMA consensus whose HA1 and HA2 frames differ at
the site, which is the outcome the rule has to catch.

| | R10, site readable | R10, frameshifted | R9 |
|---|---|---|---|
| support at the site | 0.61 | 0.42 | 0.27 |
| control support | 0.79 | 0.55 | 0.48 |

The rule above flags **all 12** of the frameshifted samples and 11 of
the 92 others, leaving 81 of the 104 firm.

Two findings decided its shape.

**The control windows carry the chemistry.** Every one of the 10 R9
samples falls below the control floor (their highest is 0.56) and the
R10 samples sit well above it (median 0.79). A rule naming R9 was
therefore doing nothing the reads were not already saying, and it is
gone. That matters beyond tidiness: a bad R10 flow cell, or a chemistry
that does not exist yet, is judged on the reads in front of it rather
than on a label.

**The ratio of the two does not work.** The obvious idea — divide the
site's agreement by the control's, so that each sample calibrates
itself — was tried and is worse: at its best it catches 8 of the 12
against 12, and at thresholds that catch more it flags most of the good
samples too. The reason is in the table above. The samples that fail are
not samples where the site is unusually hard; they are samples that are
thin, or hold two viruses, and that degrades the control just as much as
the site, so dividing one by the other cancels the signal. In the R10
samples the ratio is actually *higher* in the frameshifted ones (0.92)
than in the readable ones (0.77). The two measurements are worth more
apart than as a quotient.

The thresholds remain the numbers of one laboratory's runs. They are a
starting point: a laboratory whose runs differ should look at the
support and control support its own clean samples reach before trusting
the line between `firm` and `confirm`. Everything reported beside the
verdict — the residues, the reads behind them, both supports, and the
other sites found — does not depend on them.

## What other tools do

Surveyed 2026-11-27. Every tool below derives the site from an
assembled consensus; none of them reads it from the reads or reports
how many reads support it.

| tool | input it uses for the site | how it calls the site | read support? |
|---|---|---|---|
| [nf-flu](https://github.com/CFIA-NCFAD/nf-flu) (CFIA/PHAC) | its own consensus, annotated by VADR | residue −4 must be R or K, then counts R/K backwards, tolerating one non-basic residue: 1 basic = monobasic/LPAI, 2 = multibasic LP/HP, 3 or more = multibasic/HPAI | no |
| [Nextstrain avian-flu](https://github.com/nextstrain/avian-flu) | the HA alignment of consensus sequences | finds `GLFG`, takes the 12 bases before it and translates them; calls the motif present when R at 1 and 4 and R or K at 3 | no |
| [Nextclade](https://doi.org/10.1093/ve/veaf058) A(H5) datasets | consensus FASTA | a fixed window of HA residues 341–346; 4 or more R/K is flagged polybasic | no (only its general QC) |
| [VADR](https://github.com/ncbi/vadr) | FASTA | annotates the cleavage-site feature; makes no call | no |
| [IRMA](https://doi.org/10.1186/s12864-016-3030-6), [GenoFLU](https://doi.org/10.1016/j.virol.2023.109860), [LABEL](https://doi.org/10.1186/s12864-016-3030-6), [INSaFLU](https://doi.org/10.1186/s13073-024-01334-3), [FluMut](https://doi.org/10.1093/ve/veaf011), EPI2ME `wf-flu` | – | do not report the site at all | – |

Two of them are worth reading beside this one.

**Nextstrain's script anchors on the fusion peptide**, exactly as flucut
does, and translates the bases before it. The idea of finding the site
by what follows it is theirs as much as anyone's; what is done here is
to apply it to the reads, one at a time, and to count how many of them
agree.

**[DeGenRESOLVE](https://doi.org/10.1002/2211-5463.70371)** (2026) is
the nearest work. It takes the other route: it repairs the consensus,
adjudicating the indels a variant caller proposes by depth, by the
fraction of reads supporting them, and by a frame test over groups of
indels within 12 bases. On eight isolates it recovered the site in seven
against three to seven for the other tools it compared. Its own stated
limit is where this check begins: the result is "a dominant consensus,
not a record of within-sample variation", and an ambiguity survives only
when the minor allele is above about 40%. The two are complementary —
one mends the consensus, the other says what the reads said before there
was one.

It also measures what is at stake. Without indel calling, the cleavage
site was wrong in all five of its isolates whose site differs in length
from the reference, and the errors ran both ways: three highly
pathogenic H7 viruses were given a **monobasic** site, a low
pathogenicity H5 a polybasic one. A wrong site is a wrong pathotype, not
a wrong base.

The same point is made for Nanopore chemistry by
[Siegers et al. 2024](https://doi.org/10.1128/spectrum.01880-24), who
found the correct motif in 59.5% of R9 samples against 90.4% of R10,
with every R9 indel causing a frameshift which "could result in
incorrectly classified virulence with automated pipelines". They read
the motif from a Sanger consensus, with the pattern `(R|K)[^RKH](R|K)`,
and list minor variants as not addressed.

Reference laboratories do not use an automated caller for this: the
pathotype is established by RT-qPCR and Sanger sequencing of the site,
or by reading the consensus by hand. Not confirmed in this survey:
whether FluSurver reports the site, and whether BV-BRC annotates it.

## What it does not do

- It does not assign a pathotype. A multibasic site is the sequence
  determinant of high pathogenicity, not a diagnosis; pathotype is
  decided by the designated authorities, on more than one segment and
  usually with a test in birds.
- It does not correct the consensus. The majority over reads which make
  indels in this stretch is itself biased towards losing bases, so a
  site supported by few reads is a reason to confirm, not an answer.
- It does not tell apart the haemagglutinins of a mixed sample. It
  reports the sites it finds and their shares, and says when they
  disagree with the sample's own reference.
- It reads 8 codons. That covers the basic run of every reference flucut
  ships (the longest is 6 residues), but an insert longer than that, as
  some highly pathogenic H7 carry, is shown only in part. The site still
  reads as multibasic, so the verdict is unaffected; only the residues
  printed are cut short.
- `multibasic` here means the last residue is R or K and at least four
  of the last six are: a reading aid, not the definition used by any
  authority. Other work uses other rules, for example the pattern
  `(R|K)[^RKH](R|K)`.

## Where it is reported

- the HTML report and the dashboard: a card in each sample, and a column
  of the table of the HA in all samples;
- `summary.tsv`: `cleavage_site`, `cleavage_reads`, `cleavage_support`,
  `cleavage_multibasic`, `cleavage_verdict`;
- `flucut_summary.json`: the `cleavage_site` object of each section (see
  `docs/json_summary.md`).

It needs references (it is not reported with `-r off`) and an HA
reference which contains the fusion peptide: a reference which is a
coding sequence works, as the site is inside the ORF.

## In the code

| file | what it holds |
|---|---|
| `src/codons.h`, `.c` | the genetic code: translation, and finding a peptide in a sequence |
| `src/anchor-site.h`, `.c` | an anchor: building it from a reference, finding it in a read, reading the bases beside it. Not specific to this site |
| `src/cleavage-site.h`, `.c` | the anchors of the HA references, the motifs counted, and what the reads say |
| `src/trimmer.c` | which reads are read, per platform (`read_cleavage_site()`) |
| `src/counter.c` | the motifs of a group of input files, merged over the threads and saved for `--resume` |
