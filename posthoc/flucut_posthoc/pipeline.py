"""The run: map the reads, analyse them in parallel, write the results."""

import json
import multiprocessing
import os
import subprocess
import tempfile
import threading
from collections import Counter
from dataclasses import dataclass, field

import numpy as np
from loguru import logger

from .amplicons import CATEGORIES, Params
from .analysis import BatchResult, init_worker, process_batch
from .flucut_report import on_references, read_coverage
from .flucut_tags import (FLUCUT_CATEGORIES, has_flucut_tags,
                          load_flucut_reads)
from .irma import IrmaRun, lift_to_references, load_irma
from .mapping import read_groups, start_minimap2
from .offtarget import cluster_offtarget
from .primers import OLIGOS
from .references import (Reference, adapt_references, choose_references,
                         load_references, position_maps)
from .report import write_report
from .seqio import write_fasta
from .summary import build_summary

BATCH_READS = 2000
OFFTARGET_PER_BATCH = 100
ROW_HEADER = ("read\tunit\tunits\tread_length\treference\tstrand"
              "\tleft_position\tleft_primer\tleft_primer_errors\tleft_origin"
              "\tright_position\tright_primer\tright_primer_errors"
              "\tright_origin\ttemplate_length\tdeletions\tcategory"
              "\tidentity\n")


@dataclass
class Config:
    reads: list[str]
    outdir: str
    name: str
    references: str
    irma: list[str] = field(default_factory=list)
    flucut_reads: str | None = None
    flucut_summary: str | None = None
    flucut_group: str | None = None
    flucut_report: str | None = None
    threads: int = 8
    minimap2: str = "minimap2"
    samtools: str = "samtools"
    error_rate: float = 0.2
    terminal_window: int = 10
    min_deletion: int = 100
    command: str = ""


class Coverage:
    """Depth along each reference by category, from intervals."""

    def __init__(self, references: list[Reference], categories: list[str]):
        self.categories = {c: i for i, c in enumerate(categories)}
        self.diff = [np.zeros((len(categories), r.length + 1), np.int64)
                     for r in references]
        self.used = False

    def add(self, intervals: list[tuple[int, str, int, int]]) -> None:
        for ref, category, start, end in intervals:
            diff = self.diff[ref]
            size = diff.shape[1] - 1
            start, end = max(start, 0), min(end, size)
            if start < end:
                row = self.categories[category]
                diff[row, start] += 1
                diff[row, end] -= 1
                self.used = True

    def depth(self, ref: int) -> np.ndarray:
        return np.cumsum(self.diff[ref], axis=1)[:, :-1]


def _batches(sam, flucut: dict | None, irma: list[IrmaRun],
             slots: threading.Semaphore):
    """Batches of reads with what flucut and IRMA did with each; waits
    for a free slot, so that the mapper does not run ahead of the
    analysis."""
    batch = []
    for name, lines in read_groups(sam):
        parts = flucut.get(name) if flucut is not None else False
        used = tuple((run.label, run.used.get(name)) for run in irma)
        batch.append((name, lines, parts, used))
        if len(batch) >= BATCH_READS:
            slots.acquire()
            yield batch
            batch = []
    if batch:
        slots.acquire()
        yield batch


def _map_and_analyse(config: Config, references: list[Reference],
                     params: Params, ref_lengths: dict[str, int], lift: dict,
                     flucut: dict | None, irma: list[IrmaRun], rows_path: str,
                     trimmed: bool
                     ) -> tuple[Counter, Coverage, Coverage, list]:
    counts: Counter = Counter()
    measured = Coverage(references, CATEGORIES)
    predicted = Coverage(references, FLUCUT_CATEGORIES)
    offtarget: list = []
    workers = max(1, config.threads)
    slots = threading.Semaphore(4 * workers)
    with tempfile.TemporaryDirectory() as tmp:
        ref_path = os.path.join(tmp, "references.fasta")
        write_fasta(ref_path, [(r.name, r.seq) for r in references])
        mapper = start_minimap2(config.minimap2, ref_path, config.reads,
                                config.threads)
        gzip = subprocess.Popen(["gzip", "-c", "-3"], stdin=subprocess.PIPE,
                                stdout=open(rows_path, "wb"), text=True)
        gzip.stdin.write(ROW_HEADER)
        with multiprocessing.Pool(workers, init_worker,
                                  (references, params, ref_lengths,
                                   lift, trimmed)) as pool:
            results = pool.imap(process_batch,
                                _batches(mapper.stdout, flucut, irma, slots))
            for result in results:
                slots.release()
                counts.update(result.counts)
                measured.add(result.coverage)
                predicted.add([p for p in result.predicted
                               if p[1] in predicted.categories])
                offtarget.extend(result.offtarget[:OFFTARGET_PER_BATCH])
                gzip.stdin.write("".join(result.rows))
                if counts["reads"] % 200000 < BATCH_READS:
                    logger.info(f"{counts['reads']} reads analysed")
        gzip.stdin.close()
        gzip.wait()
        if mapper.wait() != 0:
            raise RuntimeError("minimap2 failed")
    return counts, measured, predicted, offtarget


def _load_flucut_group(path: str, group: str | None, name: str
                       ) -> dict | None:
    """The section of a group in a JSON summary of flucut: the group
    named, or the one whose name the sample name contains, or all reads
    if the run had one group."""
    with open(path) as handle:
        data = json.load(handle)
    groups = [g for g in data.get("groups", []) if g.get("processed")]
    for g in groups:
        if g["name"] == (group or name):
            return g
    if group is None:
        matching = [g for g in groups if g["name"] in name
                    or name in g["name"]]
        if len(matching) == 1:
            return matching[0]
        if len(groups) <= 1:
            return data.get("all")
    raise ValueError(f"group {group or name!r} not found in {path}; "
                     f"groups: {', '.join(g['name'] for g in groups)}")


def run(config: Config) -> dict:
    """Run the analysis of one sample; returns the summary written."""
    os.makedirs(config.outdir, exist_ok=True)
    params = Params(OLIGOS, config.error_rate,
                    config.terminal_window,
                    min_deletion=config.min_deletion)
    all_refs = load_references(config.references)
    ref_lengths = {r.name: r.length for r in all_refs}
    logger.info("choosing the reference of each segment")
    references, evidence = choose_references(
        all_refs, config.reads, config.minimap2, config.threads)
    references = adapt_references(references, config.reads,
                                  config.minimap2, config.threads)
    logger.info("references: " + ", ".join(
        f"{r.gene} {r.name}" for r in references))

    flucut, flucut_source = None, None
    if config.flucut_reads:
        flucut_source = config.flucut_reads
    elif all(has_flucut_tags(path) for path in config.reads):
        flucut_source = config.reads[0] if len(config.reads) == 1 else None
    if flucut_source:
        logger.info(f"reading the tags of flucut from {flucut_source}")
        flucut = load_flucut_reads(flucut_source,
                                   by_record=not config.flucut_reads)
    trimmed_input = flucut is not None and not config.flucut_reads
    group = None
    if config.flucut_summary:
        group = _load_flucut_group(config.flucut_summary,
                                   config.flucut_group, config.name)
    irma = []
    for path in config.irma:
        logger.info(f"reading the IRMA results in {path}")
        label = os.path.basename(path.rstrip("/")).replace(config.name, "")
        irma.append(load_irma(path, config.samtools,
                              label.strip("_") or "IRMA"))
    lifted = [lift_to_references(run_, references, config.minimap2)
              for run_ in irma]

    logger.info("mapping and analysing the reads")
    rows_path = os.path.join(config.outdir, "reads.tsv.gz")
    with_flucut = flucut is not None or config.flucut_summary or \
        config.flucut_report
    lift = position_maps(references, all_refs, config.minimap2) \
        if with_flucut else {}
    shown = None
    if config.flucut_report:
        shown = on_references(
            read_coverage(config.flucut_report,
                          config.flucut_group or config.name),
            references, lift)
    for ref in references:
        ref_lengths[ref.name] = ref.length
    counts, measured, predicted, offtarget = _map_and_analyse(
        config, references, params, ref_lengths, lift, flucut, irma,
        rows_path, trimmed_input)
    clusters = cluster_offtarget(offtarget)
    if clusters:
        write_fasta(os.path.join(config.outdir, "offtarget_amplicons.fasta"),
                    [(f"offtarget_{i} sampled_reads={c['sampled_reads']} "
                      f"template_length={c['template_length']}",
                      c["representative"])
                     for i, c in enumerate(clusters, 1)])
    summary = build_summary(
        config, params, references, evidence, counts, measured,
        predicted if predicted.used else None, clusters, irma, lifted,
        group, flucut is not None, trimmed_input, lift, shown)
    path = os.path.join(config.outdir, "posthoc_summary.json")
    with open(path + ".tmp", "w") as handle:
        json.dump(summary, handle)
    os.replace(path + ".tmp", path)
    write_report(summary, os.path.join(config.outdir, "posthoc_report.html"))
    logger.info(f"results written to {config.outdir}")
    return summary
