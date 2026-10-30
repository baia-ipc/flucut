"""Screen FASTA for the HA fusion-peptide motif GLFGA[IK]AG in 6 frames.
Report records whose motif residue is K (artefact) or other non-I.
Usage: python3 -I screen_motif.py in.fasta > out.tsv
TSV: accession, residue(I/K/other), frame, nt_start(1-based,of GLFGA..codon),
     peptide(GLFGAxAG), defline"""
import sys, re
COD = {}
B = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(AA):
  COD[B[i >> 4] + B[(i >> 2) & 3] + B[i & 3]] = a
def tr(s):
  return "".join(COD.get(s[i:i+3], "X") for i in range(0, len(s) - 2, 3))
comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
pat = re.compile(r"GLFGA(.)AG")
def recs(fn):
  acc = seq = None
  defl = ""
  for line in open(fn):
    if line.startswith(">"):
      if acc: yield acc, defl, "".join(seq)
      defl = line[1:].strip(); acc = defl.split()[0]; seq = []
    elif acc:
      seq.append(line.strip())
  if acc: yield acc, defl, "".join(seq)
print("accession\tresidue\tframe\tnt_start\tpeptide\tdefline")
for acc, defl, seq in recs(sys.argv[1]):
  seq = seq.upper()
  rc = seq.translate(comp)[::-1]
  for strand, s in (("+", seq), ("-", rc)):
    for fr in range(3):
      p = tr(s[fr:])
      for m in pat.finditer(p):
        res = m.group(1)
        aa_i = m.start()          # index in protein of 'G' of GLFGA
        nt = fr + aa_i * 3 + 1     # 1-based nt of that G on strand s
        # nt of the variable codon (6th residue): +5 residues
        print(acc, res, strand + str(fr), fr + (aa_i + 5) * 3 + 1,
              m.group(0), defl, sep="\t")
