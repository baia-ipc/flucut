# (c) 2026 Giorgio Gonnella
"""Download influenza A segment reference sequences from NCBI.

This script fetches nucleotide reference sequences for the eight
influenza A genome segments (PB2, PB1, PA, HA, NP, NA, MP, NS) from
NCBI using the public E-utilities API (esearch + efetch against
nuccore) and writes them to a FASTA file. Alongside the FASTA it
writes a companion provenance table (``<output>.provenance.tsv``)
recording, for each sequence, its accession, description, length,
the source URL it was fetched from, and the fetch date.

flucut aligns sequencing reads against these references, so having a
curated, reproducible set of references with recorded provenance
matters. By default the script prefers NCBI RefSeq influenza records
(the curated reference set) when ``--refseq-only`` is given.

The script uses the Python standard library only (no requests, no
Biopython), honoring the ``HTTPS_PROXY`` / ``https_proxy``
environment variables automatically via urllib's ProxyHandler.

Example commands:

  # A couple of RefSeq HA references for H5N1
  python3 scripts/fetch_references.py --segment HA --subtype H5N1 \\
      --max 2 --refseq-only --output references/ha_h5n1.fasta

  # All eight segments, 3 references each, into one file
  python3 scripts/fetch_references.py --segment all --max 3 \\
      --refseq-only --output references/influenza_a_all.fasta

  # With an NCBI API key (raises the rate limit)
  python3 scripts/fetch_references.py --segment NA --max 5 \\
      --api-key YOUR_KEY_HERE
"""

import argparse
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# The eight influenza A genome segments, in segment order (1-8).
SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]

# The gene names of the segments whose records NCBI annotates by the
# names of their genes rather than of the segment: segment 7 (MP)
# encodes M1 and M2, segment 8 (NS) NS1 and NEP (formerly NS2);
# MP[Gene] and NS[Gene] match none of the RefSeq records. The other
# segments are queried by their own name.
GENE_ALIASES = {
  "MP": ["M1", "M2"],
  "NS": ["NS1", "NEP", "NS2"],
}

# Be polite to NCBI: at most ~3 req/s without a key, ~10 with a key.
# We sleep a little more than the minimum to stay well within limits.
SLEEP_NO_KEY = 0.4
SLEEP_WITH_KEY = 0.15

# Network timeout in seconds for each request.
TIMEOUT = 60


def log(msg):
  """Print a progress/error message to stderr."""
  sys.stderr.write(msg + "\n")
  sys.stderr.flush()


def gene_term(segment):
  """The esearch term selecting the records of the genes of a segment.

  For example HA[Gene], or (M1[Gene] OR M2[Gene]) for MP.
  """
  genes = GENE_ALIASES.get(segment)
  if genes is None:
    return "%s[Gene]" % segment
  return "(%s)" % " OR ".join("%s[Gene]" % g for g in genes)


def build_query(segment, subtype, refseq_only):
  """Build the esearch term for one influenza A segment.

  The base query restricts to influenza A virus nucleotide records
  for the genes of the given segment (see gene_term()). Example
  resulting term:

    Influenza A virus[Organism] AND HA[Gene] AND refseq[filter]

  With a subtype (e.g. H5N1) the subtype string is added as a free
  text term to refine the search, and ``srcdb_refseq[properties]``
  (equivalent to ``refseq[filter]``) restricts to the curated RefSeq
  set when ``--refseq-only`` is requested.
  """
  parts = ['Influenza A virus[Organism]']
  parts.append(gene_term(segment))
  if subtype:
    parts.append('%s' % subtype)
  if refseq_only:
    parts.append('refseq[filter]')
  return " AND ".join(parts)


def http_get(url):
  """Perform an HTTP GET through the environment proxy, return bytes.

  urllib automatically honors HTTPS_PROXY / https_proxy via the
  default opener's ProxyHandler, which reads
  urllib.request.getproxies(). We build an opener explicitly from
  those proxies so the behavior is clear and testable.
  """
  proxies = urllib.request.getproxies()
  handler = urllib.request.ProxyHandler(proxies)
  opener = urllib.request.build_opener(handler)
  req = urllib.request.Request(
    url, headers={"User-Agent": "flucut-fetch-references/1.0"})
  with opener.open(req, timeout=TIMEOUT) as resp:
    return resp.read()


def eutils_url(tool, params):
  """Compose a full E-utilities URL for a tool with query params."""
  query = urllib.parse.urlencode(params)
  return "%s%s.fcgi?%s" % (EUTILS_BASE, tool, query)


def esearch(term, retmax, api_key):
  """Run esearch against nuccore, return a list of UID strings."""
  params = {
    "db": "nuccore",
    "term": term,
    "retmax": str(retmax),
    "retmode": "json",
  }
  if api_key:
    params["api_key"] = api_key
  url = eutils_url("esearch", params)
  try:
    raw = http_get(url)
  except urllib.error.URLError as err:
    raise RuntimeError("esearch request failed: %s" % err)
  import json
  try:
    data = json.loads(raw.decode("utf-8"))
  except ValueError as err:
    raise RuntimeError("could not parse esearch JSON: %s" % err)
  result = data.get("esearchresult", {})
  if "ERROR" in result:
    raise RuntimeError("esearch error: %s" % result["ERROR"])
  return result.get("idlist", [])


def efetch_fasta(uids, api_key):
  """Run efetch for the given UIDs, return FASTA text."""
  params = {
    "db": "nuccore",
    "id": ",".join(uids),
    "rettype": "fasta",
    "retmode": "text",
  }
  if api_key:
    params["api_key"] = api_key
  url = eutils_url("efetch", params)
  try:
    raw = http_get(url)
  except urllib.error.URLError as err:
    raise RuntimeError("efetch request failed: %s" % err)
  return raw.decode("utf-8")


def parse_fasta(text):
  """Parse FASTA text into a list of (header, sequence) tuples.

  The header is the full description line without the leading '>'.
  The sequence is the concatenated residues (newlines stripped).
  """
  records = []
  header = None
  seq_lines = []
  for line in text.splitlines():
    if line.startswith(">"):
      if header is not None:
        records.append((header, "".join(seq_lines)))
      header = line[1:].strip()
      seq_lines = []
    else:
      seq_lines.append(line.strip())
  if header is not None:
    records.append((header, "".join(seq_lines)))
  return records


def accession_from_header(header):
  """Extract the accession (first token) from a FASTA header."""
  return header.split()[0] if header else ""


def description_from_header(header):
  """Extract the description (header minus the accession token)."""
  parts = header.split(None, 1)
  return parts[1] if len(parts) > 1 else ""


def sleep_for(api_key):
  """Sleep between requests to respect NCBI rate limits."""
  time.sleep(SLEEP_WITH_KEY if api_key else SLEEP_NO_KEY)


def fetch_segment(segment, subtype, max_seqs, refseq_only, api_key):
  """Fetch references for one segment.

  Returns a list of dicts with keys: accession, description,
  sequence, source_url. Raises RuntimeError on failure.
  """
  term = build_query(segment, subtype, refseq_only)
  log("[%s] esearch: %s" % (segment, term))
  uids = esearch(term, max_seqs, api_key)
  sleep_for(api_key)
  if not uids:
    log("[%s] WARNING: no records found for this query" % segment)
    return []
  log("[%s] found %d record(s), fetching FASTA" % (segment, len(uids)))
  fasta_text = efetch_fasta(uids, api_key)
  sleep_for(api_key)
  records = parse_fasta(fasta_text)
  if not records:
    log("[%s] WARNING: efetch returned no sequences" % segment)
    return []
  efetch_params = {
    "db": "nuccore",
    "id": ",".join(uids),
    "rettype": "fasta",
    "retmode": "text",
  }
  source_url = eutils_url("efetch", efetch_params)
  out = []
  for header, seq in records:
    out.append({
      "accession": accession_from_header(header),
      "description": description_from_header(header),
      "sequence": seq,
      "source_url": source_url,
    })
  return out


def write_fasta(records, path):
  """Write records to a FASTA file, 70 residues per line."""
  with open(path, "w") as handle:
    for rec in records:
      handle.write(">%s %s\n" % (rec["accession"], rec["description"]))
      seq = rec["sequence"]
      for i in range(0, len(seq), 70):
        handle.write(seq[i:i + 70] + "\n")


def write_provenance(records, path):
  """Write the provenance TSV companion file."""
  fetch_date = datetime.now(timezone.utc).isoformat()
  header = "accession\tdescription\tlength\tsource_url\tfetch_date\n"
  with open(path, "w") as handle:
    handle.write(header)
    for rec in records:
      row = "%s\t%s\t%d\t%s\t%s\n" % (
        rec["accession"],
        rec["description"],
        len(rec["sequence"]),
        rec["source_url"],
        fetch_date,
      )
      handle.write(row)


def default_output(segment):
  """Default FASTA output path for a given segment selection."""
  name = "all" if segment == "all" else segment
  return "references/influenza_a_%s.fasta" % name


def ensure_parent_dir(path):
  """Create the parent directory of path if it does not exist."""
  import os
  parent = os.path.dirname(os.path.abspath(path))
  if parent and not os.path.isdir(parent):
    os.makedirs(parent, exist_ok=True)


def parse_args(argv):
  """Parse command-line arguments."""
  desc = "Download influenza A segment references from NCBI."
  parser = argparse.ArgumentParser(description=desc)
  seg_help = ("segment name (one of %s) or 'all'"
              % ", ".join(SEGMENTS))
  parser.add_argument(
    "--segment", default="HA",
    help=seg_help + " (default: HA)")
  parser.add_argument(
    "--subtype", default=None,
    help="optional subtype to refine the search, e.g. H5N1")
  parser.add_argument(
    "--max", dest="max_seqs", type=int, default=5,
    help="max sequences per segment (default: 5)")
  parser.add_argument(
    "--refseq-only", action="store_true",
    help="restrict to curated NCBI RefSeq influenza records")
  parser.add_argument(
    "--output", default=None,
    help="FASTA output path "
         "(default: references/influenza_a_<segment>.fasta)")
  parser.add_argument(
    "--api-key", default=None,
    help="optional NCBI API key (raises the rate limit)")
  return parser.parse_args(argv)


def resolve_segments(segment):
  """Validate and expand the --segment value to a list of names."""
  if segment == "all":
    return list(SEGMENTS)
  seg = segment.upper()
  if seg not in SEGMENTS:
    raise ValueError(
      "unknown segment %r; choose one of %s or 'all'"
      % (segment, ", ".join(SEGMENTS)))
  return [seg]


def main(argv=None):
  """Program entry point. Returns an exit code."""
  args = parse_args(argv if argv is not None else sys.argv[1:])

  try:
    segments = resolve_segments(args.segment)
  except ValueError as err:
    log("error: %s" % err)
    return 2

  if args.max_seqs < 1:
    log("error: --max must be at least 1")
    return 2

  output = args.output or default_output(args.segment)
  prov_path = output + ".provenance.tsv"

  all_records = []
  for seg in segments:
    try:
      recs = fetch_segment(
        seg, args.subtype, args.max_seqs,
        args.refseq_only, args.api_key)
    except RuntimeError as err:
      log("error while fetching segment %s: %s" % (seg, err))
      return 1
    all_records.extend(recs)

  if not all_records:
    log("error: no sequences were downloaded; nothing written")
    return 1

  try:
    ensure_parent_dir(output)
    write_fasta(all_records, output)
    write_provenance(all_records, prov_path)
  except OSError as err:
    log("error writing output files: %s" % err)
    return 1

  accs = ", ".join(r["accession"] for r in all_records)
  log("wrote %d sequence(s) to %s" % (len(all_records), output))
  log("accessions: %s" % accs)
  log("provenance written to %s" % prov_path)
  return 0


if __name__ == "__main__":
  sys.exit(main())
