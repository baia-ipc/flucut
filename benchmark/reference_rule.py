# (c) 2026 Giorgio Gonnella
"""The rule of flucut for the reference of each segment, from the votes.

A copy of auto_settings_decide_references() (src/auto-settings.c) whose
three constants can be changed, to see what other values would choose
from the votes of a run: min_votes (AUTO_REFERENCE_MIN_VOTES),
min_share (AUTO_REFERENCE_MIN_SHARE) and close_share
(AUTO_REFERENCE_CLOSE_SHARE; None: the common subtype is never
followed where a reference leads without a clear lead).
"""

DEFAULTS = {"min_votes": 20, "min_share": 0.6, "close_share": 0.5}


def leader(votes):
  """The reference with most votes (of several, the first)."""
  best = None
  for ref, n in votes.items():
    if best is None or n > votes[best]:
      best = ref
  return best


def is_clear(votes, min_votes, min_share):
  best = leader(votes)
  return (votes[best] >= min_votes and
          votes[best] >= min_share * sum(votes.values()))


def common_subtype(segments, subtype, min_votes, min_share):
  """The subtype chosen in most clear segments (then: with most votes
  in them; then: of the first segment), None if no segment is clear."""
  chosen = {}
  for s, votes in sorted(segments.items()):
    if is_clear(votes, min_votes, min_share):
      best = leader(votes)
      chosen[s] = (subtype[best], votes[best])
  common, common_key = None, (0, 0)
  for s in sorted(chosen):
    same = [v for st, v in chosen.values() if st == chosen[s][0]]
    if (len(same), sum(same)) > common_key:
      common, common_key = chosen[s][0], (len(same), sum(same))
  return common


def decide(segments, subtype, min_votes=20, min_share=0.6,
           close_share=0.5):
  """The reference chosen for each segment and how: segments maps a
  segment to the votes of its references (in the order of the reference
  set), subtype a reference to its subtype. Returns, by segment,
  (reference, how); the reference is None where none is chosen (too few
  votes: flucut then uses all); how is one of clear, follows (the
  common subtype), leader (no clear lead, not followed), none."""
  common = common_subtype(segments, subtype, min_votes, min_share)
  result = {}
  for s, votes in segments.items():
    best = leader(votes)
    if is_clear(votes, min_votes, min_share):
      result[s] = (best, "clear")
      continue
    few = votes[best] < min_votes
    follow = None
    for ref, n in votes.items():
      # (the HA and the NA, segments 4 and 6, never follow)
      if subtype[ref] != common or common is None or s in (4, 6):
        continue
      if not few and (close_share is None or
                      n < close_share * votes[best]):
        continue
      if follow is None or n > votes[follow]:
        follow = ref
    if follow is not None:
      result[s] = (follow, "follows")
    elif few:
      result[s] = (None, "none")
    else:
      result[s] = (best, "leader")
  return result
