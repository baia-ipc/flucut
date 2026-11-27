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

Two numbers decide: how many reads gave the site, and what share of them
agree on the same residues (the **support**). How much support is enough
depends on how often the reads make an indel in this stretch, which
differs between the kinds of data, so each is judged by its own number:

| data | support needed to be firm |
|---|---|
| Illumina | 0.60 |
| Nanopore R10 | 0.45 |
| Nanopore R9 | never firm on its own |

The chemistry is taken from the basecalling model in the titles of the
reads (`dna_r10.4.1_...`, `dna_r9.4.1_...`); when the titles do not say,
the R10 number is used and the report says the chemistry is unknown.

The verdict is one of

- **firm**: enough reads, enough of them agree, and they read like the
  reference of the sample;
- **confirm**: the residues reported are the majority, but too few of
  the reads agree for them to stand on their own, or the reads do not
  read like the reference. Confirm before reporting a pathotype.
- **too few reads**: fewer than 20 reads reach the site; nothing is said
  about it.

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

The two Nanopore numbers were measured; the Illumina one was not.

The Nanopore numbers come from 104 samples of avian H5N1, H5N2 and
H9N2 (94 R10, 10 R9) of one laboratory and one protocol, with the
majority taken as above:

| | R10 | R9 |
|---|---|---|
| support, median | 0.55 | 0.26 |
| support, range | 0.25 – 1.00 | 0.15 – 0.42 |

No R9 sample reached 0.45, and the R9 reads carry an indel over the site
far more often than the R10 ones, which is why R9 is never firm.

On the same samples, 12 had an IRMA consensus whose HA1 and HA2 frames
differ at the site. The verdict above was **confirm** for 10 of them and
**too few reads** for 2: none was called firm. Of the 94 R10 samples,
68 were called firm.

The Illumina number is **not** measured on a run: no Illumina data set
of this kind was at hand. It follows from the error rate instead. An
Illumina read is wrong at about one base in a thousand, so nearly every
read spells the 24 bases correctly and a clean sample should agree far
above 0.9; 0.60 is therefore deliberately permissive, and a sample which
falls below it has something else wrong with it. Treat it as a bound to
be replaced once measured, not as a calibrated threshold.

Both of these are the numbers of one laboratory's runs. They are a
starting point, not a constant of nature: a laboratory whose runs differ
should look at the support its own clean samples reach before trusting
the line between `firm` and `confirm`. What does not depend on them is
everything reported beside the verdict — the residues, the reads behind
them, the support itself and the other sites found.

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
