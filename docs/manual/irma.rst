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
- **Keeping the termini.** With ``--trim terminus`` and the elongation,
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
