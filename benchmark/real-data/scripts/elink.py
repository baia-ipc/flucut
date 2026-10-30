"""elink a list of ids between two NCBI databases, in batches (one id= per
input, so that the links are kept per input id).
Usage: python3 -I elink.py <dbfrom> <db> <linkname> <ids.txt> > out.tsv
Output: TSV fromid, toid (UIDs)."""
import sys, time, urllib.parse, urllib.request, xml.etree.ElementTree as ET
dbfrom, db, linkname, fn = sys.argv[1:5]
ids = [l.strip() for l in open(fn) if l.strip()]
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
for i in range(0, len(ids), 100):
  batch = ids[i:i+100]
  data = urllib.parse.urlencode([("dbfrom", dbfrom), ("db", db),
      ("linkname", linkname)] + [("id", x) for x in batch]).encode()
  for attempt in range(5):
    try:
      x = urllib.request.urlopen(E, data, timeout=120).read()
      break
    except Exception as e:
      print("retry", e, file=sys.stderr); time.sleep(3 * (attempt + 1))
  root = ET.fromstring(x)
  for ls in root.findall("LinkSet"):
    frm = ls.findtext("IdList/Id")
    tos = [l.findtext("Id") for l in ls.findall("LinkSetDb/Link")]
    for t in tos or ["-"]:
      print(frm, t, sep="\t")
  sys.stdout.flush()
  time.sleep(0.4)
