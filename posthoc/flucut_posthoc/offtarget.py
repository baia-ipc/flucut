"""The amplicons which are not influenza.

The reads which do not map to the references but begin with the primer
are amplicons of another template (of the host, of another organism).
Amplicons of the same site begin with the same bases after the primer:
the reads are grouped by these bases, and for each group a read of
typical length is given as its representative, to be identified by a
database search.
"""

from collections import Counter

import edlib


def cluster_offtarget(samples: list[tuple[str, int, bool, str]],
                      max_clusters: int = 10) -> list[dict]:
    """Group a sample of reads not mapped (first template bases, template
    length, both primers found, template) by their first bases; groups
    whose first bases differ by at most 3 edits are merged. Returns the
    largest groups: share of the sample, usual template length of the
    reads with both primers, and a representative template."""
    if not samples:
        return []
    keys = Counter(key for key, _, _, _ in samples)
    centres: list[str] = []
    assigned: dict[str, str] = {}
    for key, _ in keys.most_common(2000):
        for centre in centres[:50]:
            if edlib.align(key, centre, k=3)["editDistance"] >= 0:
                assigned[key] = centre
                break
        else:
            centres.append(key)
            assigned[key] = key
    members: dict[str, list[tuple[int, bool, str]]] = {}
    for key, length, both, template in samples:
        centre = assigned.get(key)
        if centre is not None:
            members.setdefault(centre, []).append((length, both, template))
    clusters = []
    ranked = sorted(members.items(), key=lambda kv: -len(kv[1]))
    for centre, reads in ranked[:max_clusters]:
        complete = sorted(r for r in reads if r[1]) or sorted(reads)
        median = complete[len(complete) // 2]
        clusters.append({
            "first_bases": centre,
            "sampled_reads": len(reads),
            "share_of_sample": len(reads) / len(samples),
            "with_both_primers": sum(1 for r in reads if r[1]),
            "template_length": median[0],
            "representative": median[2]})
    return clusters
