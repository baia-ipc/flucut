"""Accessions -> UIDs via esearch (batches). Usage: esearch_ids.py db acc.txt > tsv acc uid"""
import sys, time, urllib.parse, urllib.request, xml.etree.ElementTree as ET
db, fn = sys.argv[1:3]
accs = [l.strip() for l in open(fn) if l.strip()]
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
for i in range(0, len(accs), 200):
  b = accs[i:i+200]
  term = " OR ".join(a + "[accn]" for a in b)
  x = urllib.request.urlopen(E + "esearch.fcgi", urllib.parse.urlencode(
      {"db": db, "term": term, "retmax": 1000}).encode(), timeout=120).read()
  ids = [e.text for e in ET.fromstring(x).findall("IdList/Id")]
  # map back with esummary
  x = urllib.request.urlopen(E + "esummary.fcgi", urllib.parse.urlencode(
      {"db": db, "id": ",".join(ids)}).encode(), timeout=120).read()
  for d in ET.fromstring(x).iter("DocumentSummary"):
    print(d.findtext("Accession"), d.get("uid"), sep="\t")
  time.sleep(0.4)
