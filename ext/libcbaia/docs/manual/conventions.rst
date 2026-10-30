.. _conventions:

How the API works
-----------------

All modules follow the same conventions. Knowing them, the use of a
module can be guessed from its header; this chapter describes them from
the point of view of a program using the library. The complete coding
rules (formatting, naming, commits), which also apply to the projects
using the library, are in the developer guide (``docs/dev_manual.md``).

Objects
~~~~~~~

Most modules implement a type, whose structure is *opaque*: the header
only declares it (``typedef struct whitelist_t whitelist_t;``) and the
fields are private to the module. Objects are created by a constructor,
used through functions taking the object as first argument, and deleted
by a destructor:

.. code-block:: c

   heap_t *heap = heap_new(cmp_ints, 16);   /* constructor */
   heap_push(heap, &value);                  /* operations */
   int *smallest = heap_pop(heap);
   heap_delete(heap);                        /* destructor */

The names are ``<module>_new()``, ``<module>_<operation>()`` and
``<module>_delete()``; the destructor frees everything the object owns,
and accepts ``NULL`` (so that it can be called on objects which were
not created, e.g. in cleanup code).

A few *plain data structures* are defined in the headers, as their
fields are read by the callers: e.g. ``fastq_batch_t`` (a batch of
records), ``strbuf_t`` (a string buffer), ``hll_t`` (a HyperLogLog
sketch), ``sam_header_inputs_t`` (the read groups and programs of the
inputs), ``whitelist_match_t`` and ``complexity_t`` (results). They are
initialized and freed by ``_init()`` and ``_free()`` functions, or are
filled by the functions that return them.

Ownership of the memory
~~~~~~~~~~~~~~~~~~~~~~~

The header of each function says who owns the memory of the values it
returns:

* strings and arrays returned by functions such as ``strcat_alloc()``,
  ``path_join()``, ``path_basename()``, ``list_dir()`` or
  ``whitelist_codes()`` are
  allocated, and must be freed by the caller (with ``free()``);
* strings passed to constructors are *copied* when they are kept (e.g.
  the description of ``progress_new()``), unless documented otherwise;
* pointers returned by accessors (e.g. ``bam_reader_header()``) point
  into the object, and are valid as long as the object exists.

Errors
~~~~~~

The library distinguishes two kinds of errors:

**Errors of the input or of the environment** (a file which cannot be
opened, an invalid FASTQ record, exhausted memory, a full disk): the
program is interrupted with an error message, through the logger, and
the exit code ``EXIT_FAILURE``. The functions therefore do not return
error codes, and their results do not need to be checked: e.g.
``xfopen()`` always returns an open file, ``fastq_reader_read_batch()``
always returns valid records. The messages describe the problem in terms
of the input of the user, e.g.::

   ERROR | Opening file 'reads.fastq' failed, No such file or directory
   ERROR | Whitelist barcodes.txt, line 3: barcodes of different lengths

The exceptions are functions whose purpose is a check, which return a
boolean: e.g. ``file_exists()``, ``is_directory()``,
``reads_reader_file_complete()``, ``bam_reader_has_eof_marker()``.

**Programming errors** (a ``NULL`` pointer where an object is required,
an index out of range): the preconditions of the exported functions are
checked with ``assert()``. These checks are disabled in the optimized
builds (``-DNDEBUG``), so they never check the input.

As errors exit the program, the library is designed for **command line
programs**; it is not suitable for long-running services or for
libraries which must recover from invalid input.

The logger
~~~~~~~~~~

Functions which can fail, or which report information, take a
``Logger *logger`` argument (usually the last one). The logger is created
by the program (see :ref:`logger`), from the options of the user, and
passed to the modules. Any function accepts ``NULL`` as logger: then only
the errors are output, to the standard error. This is convenient in unit
tests and in small tools.

Threads
~~~~~~~

The modules have no global state (except constants). Different objects
can thus be used by different threads at the same time; an object is
not protected by locks, and is used by one thread at a time (or only
read, e.g. a whitelist shared by the threads after its construction).
The classes made for parallel processing are used with one object per
thread (e.g. a ``gzip_compressor_t`` for each thread), see
:ref:`recipes`.

Units and conventions of the values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Sizes are ``size_t`` and in bytes (e.g. in the allocation macros);
  counts are ``size_t``; values which can be negative or need a special
  value (e.g. ``-1`` for "unknown") are ``long`` or ``double``.
* Sequences are passed as a pointer and a length, and are not
  necessarily 0-terminated (e.g. a sequence in the buffer of a FASTQ
  batch); the functions which return 0-terminated strings say so.
* File names are often passed as a name and a suffix (``NULL`` if not
  needed), so that the callers do not need to concatenate strings:
  ``xfopen(basename, ".sorted", "w", "sorted output", logger)``.
* Booleans are ``bool`` (``stdbool.h``).
