# (c) 2026 Giorgio Gonnella
"""Find NCBI records usable as built-in references of flucut.

flucut locates the ends of the reads on complete segment references:
an end within --terminal-window bases of a reference end is a segment
terminus, the rest is internal priming. A reference must therefore
carry the untranslated regions (UTRs) at both ends of the segment,
where the Uni12/Uni13 primers bind; a record holding only the coding
sequence would turn genuine segment termini into internal priming.

For each segment and subtype asked for, this script searches NCBI
nuccore, reads the GenBank record of each candidate and keeps the
first one (or, with --per-subtype, several, spread over the hosts,
places and decades found) which:

  * is the requested segment of influenza A virus,
  * carries the requested subtype in its description,
  * begins and ends with the conserved termini of the segments of
    influenza A (the targets of the Uni12 and Uni13 primers), so that
    it covers the segment from end to end,
  * is not shorter than --min-length, has no ambiguous bases and is a
    natural isolate (see EXCLUDED and --exclude).

It writes the accessions of the records found, in the order of the
segments, ready for the list in src/builtin-references.c, and a
provenance table with the accession, description, length, the length
of the untranslated regions and the fetch date.

It uses the Python standard library only and honors the HTTPS_PROXY /
https_proxy environment variables, as scripts/fetch_references.py,
whose E-utilities helpers it uses.

Example commands:

  # one HA reference for each HA subtype
  python3 scripts/find_builtin_references.py --segment HA \\
      --subtype H1,H2,H3 --output references/ha.tsv

  # the eight segments of a subtype (e.g. to add it to the built-in set)
  python3 scripts/find_builtin_references.py --segment all \\
      --subtype H5N8 --output references/h5n8.tsv

  # several HA references of a subtype, of different lineages
  python3 scripts/find_builtin_references.py --segment HA \\
      --subtype H5 --per-subtype 4 --max 200 --output references/h5.tsv

  # the eight segments of one isolate, searched by its name
  python3 scripts/find_builtin_references.py --segment all \\
      --subtype H5N1 --strain "A/Cambodia" --from-year 2023 \\
      --output references/h5n1_cambodia.tsv
"""

import argparse
import re
import sys
import time
import urllib.error
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

sys.path.insert(0, __file__.rsplit("/", 1)[0])

from fetch_references import (  # noqa: E402
  SEGMENTS, SLEEP_NO_KEY, SLEEP_WITH_KEY, esearch, eutils_url, gene_term,
  http_get, log)

# The number of each segment in the records of NCBI ("segment 4").
SEGMENT_NUMBER = {name: i for i, name in enumerate(SEGMENTS, 1)}

# The conserved termini of every segment of influenza A, on the strand
# the records are deposited on: the 5' end is the target of the Uni12
# primer, the 3' end the reverse complement of the Uni13 primer. A
# reference of flucut must begin and end with them, so that the ends of
# the reads at a segment terminus are recognised as such; ONE_OFF
# mismatches are allowed (the fourth base of the 5' terminus is A or G).
# Records which are not natural isolates, or whose sequence is not
# trusted, recognised by their description.
EXCLUDED = ("chimeric", "synthetic", "recombinant", "vaccine", "mutant",
            "unverified", "construct", "cloning vector")

# The windows of the date of publication searched for each segment and
# subtype when several records are wanted (--per-subtype), so that the
# records of every decade can be reached (see find()). The date of
# collection is not a search field of the nucleotide database, while
# the date of publication follows it closely enough to spread the
# candidates over the lineages.
WINDOWS = [("1900", "2004"), ("2005", "2014"), ("2015", "2020"),
           ("2021", "2100")]

TERMINUS_5 = "AGCAAAAGCAGG"
TERMINUS_3 = "CCTTGTTTCTACT"
ONE_OFF = 1

# Shortest plausible complete segment, by segment (the coding sequence
# plus the UTRs; a record below this is a fragment).
MIN_LENGTH = {"PB2": 2300, "PB1": 2290, "PA": 2190, "HA": 1700,
              "NP": 1520, "NA": 1420, "MP": 990, "NS": 860}


def parse_args():
  p = argparse.ArgumentParser(
    description="Find complete segment references of influenza A.")
  p.add_argument("--segment", default="HA",
                 help="comma-separated segment names (%s) or 'all' "
                      "(default: HA)"
                 % ", ".join(SEGMENTS))
  p.add_argument("--subtype", required=True,
                 help="comma-separated subtypes, e.g. H5N8 or H1,H2,H3; "
                      "a subtype of the HA (H5) or of the NA (N8) alone "
                      "is allowed")
  p.add_argument("--strain", default="",
                 help="text the description must hold, e.g. 'A/Cambodia'")
  p.add_argument("--from-year", type=int, default=0,
                 help="earliest year of collection in the description")
  p.add_argument("--exclude", default="",
                 help="comma-separated texts a description must not "
                      "hold, besides %s" % ", ".join(EXCLUDED))
  p.add_argument("--min-length", type=int, default=0,
                 help="shortest record searched for, instead of the "
                      "length of a complete segment; a narrow range "
                      "leaves out the records of the coding sequence "
                      "alone, which are the majority for some subtypes")
  p.add_argument("--max", dest="max_seqs", type=int, default=60,
                 help="candidates looked at per segment (default: 60)")
  p.add_argument("--per-subtype", type=int, default=1,
                 help="records kept for each segment and subtype "
                      "(default: 1); several cover the lineages of a "
                      "subtype, as one record cannot: they are spread "
                      "over the hosts, places and decades found")
  p.add_argument("--output", default="",
                 help="provenance table (TSV); the accessions go to stdout")
  p.add_argument("--api-key", default="",
                 help="optional NCBI API key (raises the rate limit)")
  return p.parse_args()


def efetch_genbank(uids, api_key):
  """The GenBank records of the UIDs, as the XML of NCBI."""
  params = {"db": "nuccore", "id": ",".join(uids), "rettype": "gb",
            "retmode": "xml"}
  if api_key:
    params["api_key"] = api_key
  raw = http_get(eutils_url("efetch", params))
  return ET.fromstring(raw)


def cds_span(record):
  """First and last base of the coding sequences of the record
  (1-based), or None if it has none or only locations this script does
  not read. All of them are taken, since segment 7 holds M1 and M2 and
  segment 8 NS1 and NEP: M1 and NS1 end well before the end of the
  segment, and only the spliced M2 and NEP reach it. The ends of a
  spliced location, join(1..26,715..982), are its first and last base.
  """
  first, last = None, None
  for feature in record.iterfind("GBSeq_feature-table/GBFeature"):
    if feature.findtext("GBFeature_key") != "CDS":
      continue
    location = feature.findtext("GBFeature_location") or ""
    if "<" in location or ">" in location or "complement" in location:
      continue    # a partial end, or the other strand
    bases = [int(n) for n in re.findall(r"\d+", location)]
    if not bases:
      continue
    first = min(bases) if first is None else min(first, min(bases))
    last = max(bases) if last is None else max(last, max(bases))
  return (first, last) if first is not None else None


def collection_year(record):
  """The year of the strain name of the description, which is its last
  field (A/duck/Cambodia/X/2024(H5N1)), or 0. Two digits
  (A/swine/Belgium/1/83) are read as 19xx above 25 and as 20xx below
  it; the fields before the last are isolate numbers, not years."""
  text = record.findtext("GBSeq_definition") or ""
  found = re.search(r"A/([^()]+?)\(", text)
  if not found:
    return 0
  year = found.group(1).split("/")[-1].strip()
  if len(year) == 4 and year.isdigit() and 1900 < int(year) < 2100:
    return int(year)
  if len(year) == 2 and year.isdigit():
    return 1900 + int(year) if int(year) > 25 else 2000 + int(year)
  return 0


def subtypes_of(subtype):
  """The complete subtypes a request covers: H5N8 itself, H5 every H5N*,
  N8 every H*N8 (the records of NCBI name the complete subtype, so a
  search for "N8" alone finds nothing)."""
  if subtype[0] == "H" and "N" in subtype[1:]:
    return [subtype]
  if subtype[0] == "H":
    return [f"{subtype}N{n}" for n in range(1, 12)]
  return [f"H{h}{subtype}" for h in range(1, 17)]


def subtype_terms(subtype):
  """The texts which show the subtype in a description."""
  return [f"({s})" for s in subtypes_of(subtype)]


def subtype_query(subtype):
  """The esearch term selecting the records of the subtype."""
  found = subtypes_of(subtype)
  return found[0] if len(found) == 1 else "(%s)" % " OR ".join(found)


def check(record, segment, subtype, args):
  """The reasons why the record cannot be a reference (none: it can)."""
  out = []
  definition = record.findtext("GBSeq_definition") or ""
  length = int(record.findtext("GBSeq_length") or 0)
  sequence = (record.findtext("GBSeq_sequence") or "").upper()
  if "Influenza A virus" not in (record.findtext("GBSeq_organism") or ""):
    out.append("not influenza A virus")
  if f"segment {SEGMENT_NUMBER[segment]}" not in definition:
    out.append(f"not segment {SEGMENT_NUMBER[segment]}")
  if not any(t in definition for t in subtype_terms(subtype)):
    out.append(f"subtype {subtype} not in the description")
  if args.strain and args.strain.lower() not in definition.lower():
    out.append(f"{args.strain} not in the description")
  for word in EXCLUDED + tuple(w.strip().lower() for w in
                               args.exclude.split(",") if w.strip()):
    if word in definition.lower():
      out.append(f"{word} in the description")
  if args.from_year and collection_year(record) < args.from_year:
    out.append(f"collected before {args.from_year}")
  if length < MIN_LENGTH[segment]:
    out.append(f"{length} bases: shorter than a complete segment")
  if sequence and set(sequence) - set("ACGT"):
    out.append("ambiguous bases")
  if not sequence:
    out.append("no sequence in the record")
  else:
    if mismatches(sequence[:len(TERMINUS_5)], TERMINUS_5) > ONE_OFF:
      out.append(f"begins with {sequence[:len(TERMINUS_5)]}, not with the "
                 f"conserved terminus {TERMINUS_5}")
    if mismatches(sequence[-len(TERMINUS_3):], TERMINUS_3) > ONE_OFF:
      out.append(f"ends with {sequence[-len(TERMINUS_3):]}, not with the "
                 f"conserved terminus {TERMINUS_3}")
  return out


def mismatches(text, motif):
  """Bases of equally long <text> and <motif> which differ."""
  return sum(a != b for a, b in zip(text, motif)) + \
    abs(len(text) - len(motif))


def utr_lengths(record):
  span = cds_span(record)
  length = int(record.findtext("GBSeq_length") or 0)
  return (span[0] - 1, length - span[1]) if span else (0, 0)


def strain_origin(definition):
  """The host and the place of the strain name of a description
  (A/duck/Cambodia/X/2022(H5N1) gives duck and Cambodia; a human strain
  has no host field: A/Cambodia/X/2023). Empty where the name does not
  read that way."""
  found = re.search(r"A/([^()]+?)\(H\d+N\d+\)", definition)
  if not found:
    return ("", "")
  parts = [p.strip().lower() for p in found.group(1).split("/")]
  if len(parts) < 3:
    return ("", "")
  return ("human", parts[0]) if len(parts) == 3 else (parts[0], parts[1])


def spread(groups, wanted):
  """Up to <wanted> records of the candidate groups (one group per
  window of publication dates), as different from each other as the
  descriptions tell: one record of each group in turn, and within a
  group one of a host and place not taken yet."""
  chosen, seen = [], set()
  while len(chosen) < wanted:
    before = len(chosen)
    for group in groups:
      if len(chosen) >= wanted:
        break
      pick = next((r for r in group if r not in chosen and
                   strain_origin(r[1]) not in seen), None)
      if pick is None:
        pick = next((r for r in group if r not in chosen), None)
      if pick is not None:
        chosen.append(pick)
        seen.add(strain_origin(pick[1]))
    if len(chosen) == before:   # every group is exhausted
      break
  return chosen


def candidates(segment, subtype, args, window, wanted):
  """The records of the segment and subtype collected in one time
  window which can be references, and why the others cannot."""
  # the length range leaves out most records which hold the coding
  # sequence alone, without the UTRs (the check of the record decides)
  shortest = args.min_length or MIN_LENGTH[segment]
  terms = ["Influenza A virus[Organism]", gene_term(segment),
           subtype_query(subtype), "complete cds[Title]",
           f"{shortest}:{shortest + 400}[SLEN]"]
  if window:
    terms.append(f'("{window[0]}"[PDAT] : "{window[1]}"[PDAT])')
  if args.strain:
    terms.append(f'"{args.strain}"')
  uids = esearch(" AND ".join(terms), args.max_seqs, args.api_key)
  if not uids:
    return [], {}
  time.sleep(SLEEP_WITH_KEY if args.api_key else SLEEP_NO_KEY)
  refused, found = {}, []
  for batch in [uids[i:i + 20] for i in range(0, len(uids), 20)]:
    for record in efetch_genbank(batch, args.api_key):
      accession = record.findtext("GBSeq_primary-accession") or ""
      version = record.findtext("GBSeq_accession-version") or accession
      problems = check(record, segment, subtype, args)
      if problems:
        refused[version] = "; ".join(problems)
        continue
      found.append((version, record.findtext("GBSeq_definition") or "",
                    int(record.findtext("GBSeq_length") or 0),
                    utr_lengths(record), collection_year(record)))
    if len(found) >= 4 * wanted:
      break
    time.sleep(SLEEP_WITH_KEY if args.api_key else SLEEP_NO_KEY)
  return found, refused


def find(segment, subtype, args, wanted=1):
  """Up to <wanted> records of the segment and subtype which can be
  references, each as (accession, description, length, UTRs, year); the
  records kept are spread over the hosts, places and decades found.

  For several records the search is repeated in each window of WINDOWS:
  the records of a subtype published last are often one submission of
  one place and host, which would leave the lineages of the other
  decades out."""
  windows = WINDOWS if wanted > 1 else [None]
  groups, refused = [], {}
  for window in windows:
    part, why = candidates(segment, subtype, args, window, wanted)
    groups.append(part)
    refused.update(why)
  found = [r for group in groups for r in group]
  if not found:
    log(f"{segment} {subtype}: no candidate is a complete segment"
        + ("; e.g. " + "; ".join(f"{k}: {v}"
                                 for k, v in list(refused.items())[:3])
           if refused else ""))
  return spread(groups, wanted)


def main():
  args = parse_args()
  segments = SEGMENTS if args.segment == "all" else \
    [s.strip() for s in args.segment.split(",")]
  for segment in segments:
    if segment not in SEGMENTS:
      sys.exit(f"unknown segment: {segment}")
  fetched = datetime.now(timezone.utc).strftime("%Y-%m-%d")
  rows = []
  for subtype in args.subtype.split(","):
    for segment in segments:
      for found in find(segment, subtype.strip(), args, args.per_subtype):
        accession, definition, length, (utr5, utr3), year = found
        host, place = strain_origin(definition)
        log(f"{segment} {subtype}: {accession} ({length} bases, UTRs "
            f"{utr5} and {utr3}, {host or '?'} {place or '?'} "
            f"{year or '?'})")
        rows.append((subtype.strip(), segment, accession, definition,
                     length, utr5, utr3, host, place, year))
  for subtype in args.subtype.split(","):
    picked = [r[2] for r in rows if r[0] == subtype.strip()]
    print(f"/* {subtype.strip()} */")
    print(",".join(picked))
  if args.output:
    with open(args.output, "w") as handle:
      handle.write("subtype\tsegment\taccession\tdescription\tlength\t"
                   "utr5\tutr3\thost\tplace\tyear\tfetched\n")
      for row in rows:
        handle.write("\t".join(str(v) for v in row) + f"\t{fetched}\n")
    log(f"provenance written to {args.output}")
  missing = [(s.strip(), g) for s in args.subtype.split(",")
             for g in segments
             if not any(r[0] == s.strip() and r[1] == g for r in rows)]
  if missing:
    sys.exit("not found: " + ", ".join(f"{g} {s}" for s, g in missing))


if __name__ == "__main__":
  main()
