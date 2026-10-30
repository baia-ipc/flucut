"""Merge SRA runinfo CSV files (no header) into one TSV with a header.
Usage: python3 -I runinfo.py out.tsv in1.csv ...
"""
import csv, sys
COLS = ("Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,"
  "size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,"
  "LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,"
  "Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,"
  "SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,"
  "g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,"
  "Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,"
  "dbgap_study_accession,Consent,RunHash,ReadHash").split(",")
seen = set()
with open(sys.argv[1], "w") as out:
  out.write("\t".join(COLS) + "\n")
  for fn in sys.argv[2:]:
    for row in csv.reader(open(fn)):
      if not row or row[0] == "Run" or row[0] in seen:
        continue
      seen.add(row[0])
      row = (row + [""] * len(COLS))[:len(COLS)]
      out.write("\t".join(x.replace("\t", " ") for x in row) + "\n")
print(len(seen), "runs", file=sys.stderr)
