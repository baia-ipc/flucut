# (c) 2026 Giorgio Gonnella
"""Read the reference choice from the log of flucut (-v or -D).

For each segment with several references, flucut reports the reference
chosen and how (-v), the votes of each reference (-D), and the
agreement of the probes within a read (-v; see src/analysis/read-agreement.h).
"""

import re

CHOICE = re.compile(
    r"reference of segment (\d) .*?: (\S+) \((\d+) of (\d+) probes of the "
    r"first (\d+) reads(?:, \d+ for its subtype)?; (\d+) more match several "
    r"references\)(?:; (.*))?$")
NOT_CHOSEN = re.compile(
    r"reference of segment (\d) .*?: not chosen, too few probes \((\d+)\) "
    r"in the first \d+ reads \((\d+) more match several references\)")
VOTES = re.compile(r"votes of segment (\d) .*: (\S+) (\d+)$")
AGREEMENT = re.compile(
    r"probes within a read of segment (\d) .*: (\d+) reads with several "
    r"voting probes, all for one reference in (\d+) .*expected: (.*)$")
BY_REFERENCE = re.compile(r"(\S+) (\d+) \(([0-9.]+)\)")


def parse(text):
  """The choice of each segment, as a dictionary by segment number of
  dictionaries with: chosen (None if not chosen), how ("" if by its own
  votes), total (votes of the segment), ambiguous, votes (by
  reference), reads (with several voting probes), agreeing and expected
  (reads agreeing on each reference, observed and expected)."""
  segments = {}

  def entry(s):
    return segments.setdefault(int(s), {
        "chosen": None, "how": "", "total": 0, "ambiguous": 0,
        "votes": {}, "reads": 0, "agreeing": {}, "expected": {}})

  for line in text.splitlines():
    m = VOTES.search(line)
    if m:
      entry(m.group(1))["votes"][m.group(2)] = int(m.group(3))
      continue
    m = CHOICE.search(line)
    if m:
      e = entry(m.group(1))
      e["chosen"] = m.group(2)
      e["total"] = int(m.group(4))
      e["ambiguous"] = int(m.group(6))
      e["how"] = m.group(7) or ""
      continue
    m = NOT_CHOSEN.search(line)
    if m:
      e = entry(m.group(1))
      e["total"] = int(m.group(2))
      e["ambiguous"] = int(m.group(3))
      continue
    m = AGREEMENT.search(line)
    if m:
      e = entry(m.group(1))
      e["reads"] = int(m.group(2))
      for ref, observed, expected in BY_REFERENCE.findall(m.group(4)):
        e["agreeing"][ref] = int(observed)
        e["expected"][ref] = float(expected)
  return segments
