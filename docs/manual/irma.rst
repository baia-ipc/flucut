Using the reads with IRMA
=========================

The trimmed reads of each sample are the input of IRMA, e.g. with the
configuration for Nanopore reads:

::

   IRMA FLU-minion trimmed/barcode01.fastq.gz sample01

What was checked
----------------

These points were measured with IRMA 1.4.0 (FLU-minion) on simulated
reads (``benchmark/RESULTS.md``):

- **The tags.** IRMA accepts the tags in the titles: the consensus was
  the same with and without them. ``--no-tags`` writes the titles
  unchanged, if another program needs it.
- **The ends of the consensus.** With FLU-minion as it is, IRMA does not
  extend the consensus beyond its own references: 20 to 58 bases are
  missing at each end of the segments, with or without trimming. With
  the reference elongation (an external configuration file with
  ``SKIP_E=0``, given with ``-c``), the consensus of the reads trimmed
  by flucut ends exactly where the primers began: 12 bases short at one
  end and 13 at the other, without errors. With untrimmed reads, the
  consensus extends into the flanks, with insertions.
- **Keeping the termini.** With ``--references`` (which keeps the
  conserved termini at the terminal ends) and the elongation,
  the consensus has complete ends, but the bases there are the
  primer's: the fourth base is an ambiguity (R) in all segments, and
  where a segment differs from the primer, the primer's base is kept.
- **Low coverage.** FLU-minion replaces with N the bases supported by
  fewer than 50 reads (``MIN_CONS_SUPPORT``) in the amended consensus:
  with fewer reads, the amended consensus is all N, whatever the
  trimming; the majority consensus (``A_<gene>.fasta``) is still
  written.

The configuration file for the elongation:

::

   echo 'SKIP_E=0' > elong.sh
   IRMA FLU-minion trimmed/barcode01.fastq.gz sample01 -c elong.sh

Trimming against IRMA's own references
--------------------------------------

flucut can classify the read ends against the sequences IRMA assembles
with, so that both tools work on the same references. IRMA keeps one
seed sequence per segment and subtype in
``IRMA_RES/modules/<module>/reference/consensus.fasta``, but each seed
holds the coding sequence only, without the untranslated regions where
the Uni12 and Uni13 primers bind. flucut recognises a read end at a
segment terminus by its distance to the end of the reference, so a seed
used as it is would turn genuine segment termini into internal priming.

``scripts/irma_references.py`` builds references from the seeds: to each
seed it attaches the untranslated regions of a complete record of the
same segment and, for the HA and the NA, of the same subtype (the
built-in reference set, by default). The coding sequence is then IRMA's
and the ends reach the termini.

::

   python3 scripts/irma_references.py \
       --irma /opt/flu-amd/IRMA_RES --save-config

This writes ``irma-flu-references.fasta`` next to the cached built-in set
and, with ``--save-config``, stores in the configuration file of flucut
(``$XDG_CONFIG_HOME/flucut/config``, else ``~/.config/flucut/config``)
where IRMA is and which set to use::

   irma = /opt/flu-amd/IRMA_RES
   irma_module = FLU
   references = /home/you/.cache/flucut/irma-flu-references.fasta

A run without ``-r`` then uses that set instead of the built-in one;
``-r <file>`` and ``-r off`` still override it. The surgery is done once
after installing; when IRMA updates its seeds, run

::

   python3 scripts/irma_references.py --from-config

which takes the stored paths and writes the set again. ``--donors
<file>`` gives other complete records for the untranslated regions, and
``--output`` writes the set elsewhere (``--provenance`` the table of
what was joined: the seed, its subtype, the record whose untranslated
regions were taken and their lengths).

The untranslated regions of a segment are nearly the same in every
influenza A virus, but their length differs between the subtypes of the
HA and of the NA (the 5' untranslated region of the HA is 17 bases in
H14 and 43 in H2), which is why the donor of an HA or NA seed must be of
its subtype: a donor of another subtype would place the start of the
coding sequence at the wrong distance from the terminus. The seeds of
influenza B are left out, as the primer scheme is that of influenza A.

On simulated reads of 17 genomes outside the built-in set (subtypes H1
to H3, H5, H7, H9 and H10 with six NA subtypes, hosts from poultry to
horses, 1934 to 2017; ``benchmark/RESULTS.md``), this set named the
subtype of all 51 samples correctly, where the built-in set named 48,
and classified the read ends as well as it (the reads taken for
internally primed outside the HA, which the simulation does not prime,
were 0.6% against 0.5%). Its seeds cover the subtypes IRMA knows; the
built-in set covers every HA and NA subtype, which matters for a sample
of a subtype IRMA has no seed for.

.. warning::

   A reference built this way is not a record of a database: its coding
   sequence and its untranslated regions come from two different
   viruses, which its description and the provenance table state. Use it
   to trim the reads, never as a reference for an assembly.
