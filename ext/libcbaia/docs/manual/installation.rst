.. _installation:

Building and using the library
------------------------------

Requirements
~~~~~~~~~~~~

libcbaia is written in C11 with GNU extensions (``-std=gnu11``), for
POSIX systems (Linux, macOS). To build it, only a C compiler (gcc or
clang) and GNU ``make`` are needed: zlib is included in the repository
(``ext/zlib``), as is khash (``ext/khash``, the hash tables of klib,
used by some modules and available to the projects). valgrind is used by
the tests, if installed.

Building and testing
~~~~~~~~~~~~~~~~~~~~

::

   git clone https://github.com/baia-ipc/libcbaia
   cd libcbaia
   make          # lib/libcbaia.a and lib/libcbaia-zlib.a
   make test     # unit tests (under valgrind) and example programs

The targets of the Makefile:

=================  ==========================================================
``make``           build the static libraries ``lib/libcbaia.a`` and
                   ``lib/libcbaia-zlib.a`` (the included zlib)
``make test``      build and run the unit tests (``bin/cbaia-tests``,
                   under valgrind if available) and the example programs
                   of this manual
``make fuzz``      build the fuzzing harnesses (clang; see
                   :ref:`testing`)
``make install``   install the libraries and the headers (see below)
``make manual``    build this manual (``manual/libcbaia-manual.pdf``;
                   requires Sphinx and LaTeX)
``make clean``     remove the object and dependency files
``make cleanup``   also remove the libraries and the test program
=================  ==========================================================

**Compilation flags.** The flags are selected with the ``DEBUGFLAGS``
variable; predefined sets are:

==========================  ============================================
``DEBUGFLAGS_NDEBUG``       ``-O3 -DNDEBUG`` (default)
``DEBUGFLAGS_DEBUG``        ``-gdwarf-3 -g3 -O0``: debugging, valgrind
``DEBUGFLAGS_PROF``         ``-no-pie -pg`` and the debug flags: gprof
``DEBUGFLAGS_SANITIZE``     ``-fsanitize=address,undefined`` and the debug
                            flags
==========================  ============================================

e.g. ``make test DEBUGFLAGS='$(DEBUGFLAGS_SANITIZE)'``. When the flags
change, the objects are rebuilt (the flags of the last build are saved
in ``obj/.flags``).
Note that ``-DNDEBUG`` disables the ``assert()`` checks, which are used
for programming errors only, never for checking the input.

**System zlib.** ``make SYSTEM_LIBS=1`` uses the zlib of the system
instead of the included copy (e.g. to profit from its security updates,
or when packaging the library for a distribution); the zlib header files
are then required (e.g. the ``zlib1g-dev`` package of Debian/Ubuntu).

Installing the library
~~~~~~~~~~~~~~~~~~~~~~

::

   make install PREFIX=/usr/local

installs ``libcbaia.a`` and ``libcbaia-zlib.a`` in ``$PREFIX/lib``, and
the headers (with ``khash.h``) in ``$PREFIX/include/cbaia``. ``DESTDIR``
is supported for staged installations; ``make uninstall`` removes the
installed files.

A program then includes the headers with the directory, and links the
libraries::

   #include "cbaia/logger.h"
   #include "cbaia/fastq-reader.h"

::

   cc -std=gnu11 -O2 -o prog prog.c -lcbaia -lcbaia-zlib -lm

(with ``-lcbaia -lz -lm``, if the library was built with
``SYSTEM_LIBS=1``; add ``-lpthread`` for multithreaded programs).
Alternatively, add ``-I$PREFIX/include/cbaia`` and include the headers
by their name only, as in the vendored setup below.

Vendoring the library
~~~~~~~~~~~~~~~~~~~~~

The recommended way of using the library in a project is to include a
copy of it in the project repository (*vendoring*), e.g. in
``ext/libcbaia``: the users of the project then do not need to install
the library, and each project uses a known version of it.

**Adding and updating the copy.** With ``git subtree``, the copy is part
of the project repository, and can be updated from the library
repository::

   # add the library to the project
   git subtree add --prefix ext/libcbaia \
     https://github.com/baia-ipc/libcbaia main --squash

   # update the copy to the latest version of the library
   git subtree pull --prefix ext/libcbaia \
     https://github.com/baia-ipc/libcbaia main --squash

Changes to the library are made in the library repository (with its
tests), pushed, and then pulled into the projects, not made in the
copies: otherwise the copies diverge.

**The Makefile fragment.** ``cbaia.mk`` provides the variables and rules
needed by the Makefile of the project::

   CBAIA_DIR = ext/libcbaia
   include $(CBAIA_DIR)/cbaia.mk

   CFLAGS += -Wall -Wextra -std=gnu11 $(CBAIA_INCLUDES)
   LDLIBS += $(CBAIA_LDLIBS) -lpthread

   bin/prog: src/prog.c $(CBAIA_LIBS)
   	$(CC) $(CFLAGS) $(DEBUGFLAGS) -o $@ src/prog.c $(CBAIA_LIB) $(LDLIBS)

The variables:

=====================  ======================================================
``CBAIA_DIR``          the directory of the copy (set by the project)
``CBAIA_INCLUDES``     the ``-I`` options for the headers of the
                       components, khash and (unless ``SYSTEM_LIBS=1``) zlib
``CBAIA_LIB``          ``libcbaia.a``, to be linked
``CBAIA_LDLIBS``       the libraries to link after it: the included zlib
                       (or ``-lz``) and ``-lm``
``CBAIA_LIBS``         all libraries built by libcbaia, as prerequisites
                       of the programs
``CBAIA_EMBED``        the script ``dev/embed.sh``, which writes a data
                       file (e.g. an HTML template) as a C string literal,
                       to be included by a source (see below)
=====================  ======================================================

The library is built by its own Makefile, with the ``DEBUGFLAGS`` of
the project, if defined, so that e.g. a debug build of the project also
uses a debug build of the library. The Makefile of the library is run at
each build of the project, and rebuilds the library when one of its
sources, or the compilation flags, changed since the previous build (the
flags are saved in ``obj/.flags`` of the library); otherwise, the
library is left unchanged and the programs are not relinked. ``make
cbaia-clean`` removes the build of the library.

.. note::

   As the library is always built with the flags of the project, a
   project compiled with link-time optimization (``-flto``, as Lumipore)
   can inline the functions of the library, e.g. in the inner loops, also
   after the library was built with other flags (e.g. by ``make test`` in
   its directory): the next build of the project rebuilds it with
   ``-flto``.

**Linking only what is used.** As ``libcbaia.a`` is a static library,
only the modules used by a program are linked into it. In a project with
several programs and modules, the script ``dev/objdeps.sh`` computes the
object files of the *project* to link into each program, from the
headers which it includes (directly or indirectly)::

   dev/objdeps.sh src obj src/prog.c   # e.g. obj/samples-info.o ...

so that adding a module to a project only requires creating its ``.h``
and ``.c`` files and including the header.

**Data files.** Text used by the code but better kept in its own
format, e.g. an HTML template, can be embedded at compile time:
``$(CBAIA_EMBED) <file>`` writes the file as a C string literal (without
the lines made only of a comment, ``/* ... */`` or ``<!-- ... -->``)::

   obj/%-html.inc: src/%.html $(CBAIA_EMBED)
   	$(CBAIA_EMBED) $< > $@

   /* in the source (compiled with -Iobj) */
   static const char *const page_tmpl =
   #include "page-html.inc"
     ;

The library embeds its own data files (e.g. the style of the HTML pages,
``src/io/html-page.css``) in the same way.
