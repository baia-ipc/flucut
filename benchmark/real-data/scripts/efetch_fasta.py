"""Fetch FASTA for a list of UIDs via epost+efetch in batches.
Usage: python3 -I efetch_fasta.py db uids.txt out.fasta"""
import sys, time, urllib.parse, urllib.request
db, fn, out = sys.argv[1:4]
uids = [l.strip() for l in open(fn) if l.strip()]
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
with open(out, "w") as o:
  for i in range(0, len(uids), 300):
    b = uids[i:i+300]
    # epost
    for attempt in range(5):
      try:
        x = urllib.request.urlopen(E + "epost.fcgi", urllib.parse.urlencode(
            {"db": db, "id": ",".join(b)}).encode(), timeout=120).read().decode()
        import re
        we = re.search(r"<WebEnv>(.+?)</WebEnv>", x).group(1)
        qk = re.search(r"<QueryKey>(.+?)</QueryKey>", x).group(1)
        f = urllib.request.urlopen(E + "efetch.fcgi?" + urllib.parse.urlencode(
            {"db": db, "query_key": qk, "WebEnv": we, "rettype": "fasta",
             "retmode": "text"}), timeout=300).read().decode()
        o.write(f); break
      except Exception as e:
        print("retry", i, e, file=sys.stderr); time.sleep(3 * (attempt + 1))
    print("done", i + len(b), file=sys.stderr); o.flush()
    time.sleep(0.4)
