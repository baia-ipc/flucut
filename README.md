# libcbaia

libcbaia is a library of generic C modules for command line programs, in
particular bioinformatics tools working with sequencing reads. It is shared
by the C projects of BAIA at the Institut Pasteur du Cambodge (e.g.
[Lumipore](https://github.com/baia-ipc/lumipore)). The modules are small,
follow common conventions and depend only on the C standard library, POSIX
and zlib.

The modules are organized in four components:

- `core`: memory allocation, logging, command line interface, files,
  strings, numeric utilities, progress display, unit tests;
- `ds`: data structures (Bloom filter, set, heap, HyperLogLog);
- `io`: external sorting, parallel gzip and BGZF compression, MD5, HTML
  reports;
- `seq`: FASTQ and BAM reading and writing, SAM tags and headers, edit
  distance of short sequences, k-mer index of sequences, barcode
  whitelists, cell calling, library complexity.

## Installation

Only a C compiler (C11 with GNU extensions) and `make` are required: zlib
and khash are included in the repository (`ext`).

```
git clone https://github.com/baia-ipc/libcbaia
cd libcbaia
make
```

This creates the static libraries `lib/libcbaia.a` and
`lib/libcbaia-zlib.a`; `make SYSTEM_LIBS=1` uses the zlib of the system
instead. The tests (unit tests under valgrind, and the example programs of
the manual) are run with `make test`.

## Usage

Install the library and its headers (`make install PREFIX=/usr/local`), then
include e.g. `"cbaia/logger.h"` and link with
`-lcbaia -lcbaia-zlib -lm` (`-lcbaia -lz -lm` with `SYSTEM_LIBS=1`).

Alternatively, include a copy of the repository in the project (e.g. in
`ext/libcbaia`, with `git subtree`) and use the Makefile fragment `cbaia.mk`,
which builds the library as needed:

```
CBAIA_DIR = ext/libcbaia
include $(CBAIA_DIR)/cbaia.mk

CFLAGS += $(CBAIA_INCLUDES)
LDLIBS += $(CBAIA_LDLIBS)

bin/prog: src/prog.c $(CBAIA_LIBS)
	$(CC) $(CFLAGS) -o $@ src/prog.c $(CBAIA_LIB) $(LDLIBS)
```

The modules are then included by their name (e.g. `#include "logger.h"`).
Complete example programs (a read counter, a barcode demultiplexer, a
logging demo, an HTML report) are in [docs/manual/examples](docs/manual/examples).

## Documentation

The user manual, for the developers using the library or contributing to it,
is [manual/libcbaia-manual.pdf](manual/libcbaia-manual.pdf): building and
including the library, a first program, the conventions of the API, each
module with examples, recipes, the tests, and the API reference. Its source
is in [docs/manual](docs/manual) (Sphinx); `make manual` rebuilds the PDF
(requires Sphinx and LaTeX). The conventions of the code are described in
the developer guide ([docs/dev_manual.md](docs/dev_manual.md)).

## License

libcbaia is released under the ISC license (see the `LICENSE` file).
