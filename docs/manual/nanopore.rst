.. _nanopore:

Nanopore
--------

.. _minknow-metadata:

The ``nanopore`` directory contains modules specific to Oxford Nanopore
sequencing data. ``minknow-metadata`` converts metadata fields in MinKNOW
read titles to SAM tags. Its API and field mapping are described in the
:ref:`read metadata section <minknow-metadata>` of the sequences chapter.

Read metadata of MinKNOW (``minknow-metadata.h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The FASTQ files written by MinKNOW contain the read metadata in the
title as fields ``key=value``, separated by spaces (e.g. ``runid=...
read=12 ch=7 start_time=... flow_cell_id=... barcode=barcode01``), while
Dorado writes SAM tags. ``minknow_metadata_append_sam_tag(&out, field,
len)`` converts a field to a SAM tag (appended with a tab; returns
``false``, without changing the output, for other fields):

.. list-table::
   :header-rows: 1
   :widths: 30 12 58

   * - MinKNOW field
     - SAM tag
     - notes
   * - ``runid``
     - ``ri:Z``
     - 
   * - ``read``
     - ``rn:i``
     - as Dorado
   * - ``ch``
     - ``ch:i``
     - as Dorado
   * - ``start_time``
     - ``st:Z``
     - as Dorado
   * - ``flow_cell_id``
     - ``PU:Z``
     - standard tag (platform unit), as in the read groups of Dorado
   * - ``protocol_group_id``
     - ``pg:Z``
     - 
   * - ``sample_id``
     - ``si:Z``
     - 
   * - ``barcode``
     - ``ob:Z``
     - ``MINKNOW_BARCODE_SAM_TAG`` (ONT barcode; Dorado writes it in ``BC``)
   * - ``barcode_alias``
     - ``oa:Z``
     - 
   * - ``parent_read_id``
     - ``pi:Z``
     - as Dorado
   * - ``basecall_model_version_id``
     - ``bm:Z``
     - 

The tags which are not standard contain a lowercase letter (reserved for
end users by the SAM specification) and are not used by Dorado and
minimap2. The integer fields are written as ``Z`` if their value is not
an integer.

