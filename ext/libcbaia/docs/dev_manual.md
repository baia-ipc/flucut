# Developer guide

This guide describes how the C code of libcbaia, and of the projects
using it, is organized and written. New code shall follow these
conventions, so that the code bases stay uniform and easy to read.

## Repository layout

The library:

```
Makefile          build, test and installation targets
cbaia.mk          Makefile fragment for the projects vendoring the library
src/<component>/  sources and headers of the modules of each component:
  core/           memory allocation, logging, files, strings, numeric
                  utilities, command line, unit test macros
  ds/             data structures
  io/             external sorting, gzip and BGZF compression, MD5,
                  HTML reports
  seq/            sequences, FASTQ and BAM files, SAM tags and headers,
                  barcodes, cells, library complexity
tests/            unit test program
dev/, scripts/    scripts for the build and the projects (object lists,
                  embedded data files, program tests)
docs/             documentation (this file)
```

A project using the library:

```
Makefile          build and test targets
src/              all C sources and headers of the project (flat)
  <module>.h      interface of a module
  <module>.c      implementation of a module
  <program>.c     one file for each executable (contains main())
  filesfx.h       file suffix constants
ext/              external code (e.g. a copy of libcbaia, other libraries)
dev/              developer scripts used by the build system
scripts/          helper scripts (e.g. the test library)
testdata/         input files and expected results of the tests
docs/             documentation
bin/, obj/, .d/   build outputs (created by make, not under version control)
```

The header names are unique over all components of the library and the
projects, so that the modules include each other by name (`#include
"logger.h"`), with the directories in the include path.

File names are lower case, words separated by hyphens
(e.g. `read-layout.c`, `samples-info.h`). The name of a module is the
file name without suffix; the corresponding C identifiers use
underscores instead of hyphens (`read_layout_...`).

## Building

The build uses a single `Makefile` and the C compiler (`cc`):

- `make` builds the library (`lib/libcbaia.a`) or, in the projects, all
  programs into `bin/`; objects go into `obj/`;
- `make test` builds and runs all tests;
- `make clean` removes objects and dependency files, `make cleanup` also
  the binaries.

Header dependencies are tracked automatically (`-MMD -MP`, dependency files
in `.d/`).

Data used by the code but better written in its own format (e.g. the CSS
style and the HTML templates of the reports) is kept in data files next
to the sources and embedded at compile time: `dev/embed.sh` writes a file
as a C string literal (removing the lines made only of a comment; a CSS
file is also compacted, one rule per line, so that it can be formatted
for reading with empty lines and indentation), which the source
includes, e.g.

```c
static const char *const page_css =
#include "html-page-css.inc"
  ;
```

In the library, the files of the components with an extension listed in
`DATA_EXTS` of the Makefile (`html-page.css` gives `obj/html-page-css.inc`)
are converted automatically; a project adds a rule using `$(CBAIA_EMBED)`
(see `cbaia.mk`). The Makefiles set `.DELETE_ON_ERROR:`, so that a failed
conversion does not leave a truncated `.inc`, which the next run would
consider up to date (`cbaia.mk` does not set it, as it applies to the
whole Makefile including it: a project using `$(CBAIA_EMBED)` sets it).

The list of the components (subdirectories of `src`) is
`CBAIA_COMPONENTS` of `cbaia.mk`, also read by the Makefile of the
library and by `tests/fuzz/build.sh`; a new component is added there.

In the projects, the object files to link into each program are computed
from the headers it includes (by `dev/objdeps.sh`: each `x.h` included
directly or indirectly, also through the `.c` files of the modules, pulls
in `x.o`, if `x.c` exists). Thus, adding a module only requires creating
the `.h`/`.c` pair and including the header; programs are listed in
`PROGS`. The library is used either installed (`-lcbaia`) or vendored
(a copy in `ext/libcbaia`, built through `cbaia.mk`).

The compilation flags are selected by the `DEBUGFLAGS` variable:

| variable | flags | use |
|---|---|---|
| `DEBUGFLAGS_NDEBUG` | `-O3 -DNDEBUG -flto` | default, production and benchmarks |
| `DEBUGFLAGS_DEBUG` | `-gdwarf-3 -g3 -O0` | debugging, valgrind |
| `DEBUGFLAGS_PROF` | `-no-pie -pg` + debug flags | profiling with gprof |

e.g. `make DEBUGFLAGS='$(DEBUGFLAGS_DEBUG)'`. Note that `-DNDEBUG` disables
the `assert()` checks: they are thus used for programming errors only, never
for checking the input.

External libraries are kept to a minimum. Small single-header libraries
(e.g. for hashing) can be copied into `ext/` unchanged, with their license;
system libraries (e.g. zlib) are linked via `LDLIBS`.

## Source files

### Header of each file

Each file starts with a copyright line and, for modules, a short
description of the purpose (and, if useful, of the implementation):

```c
/* (c) 2026 Giorgio Gonnella */

/*
   Fuzzy matching of the read indices to the indices of the samples.

   Implementation:
   - exact matches are looked up in a hash table
   - otherwise, the edit distance to all indices is computed
*/
```

Headers have an include guard, named after the file in upper case:

```c
#ifndef SAMPLES_INFO_H
#define SAMPLES_INFO_H
...
#endif
```

### Formatting

- indentation by 2 spaces, no tabs, in C as in the Python and shell
  scripts (e.g. `docs/manual/*.py`, `tests/fuzz/*.sh`);
- lines up to 80 characters, in C, Python and shell, whenever possible
  (long lines are wrapped: in Python inside parentheses or with implicit
  concatenation of strings, in shell with `\`); continuation lines of
  function signatures and calls are indented by 4 spaces or aligned to
  the opening parenthesis;
- the opening brace is on the same line of the statement or function
  signature; exception: the `parseargs()` function of programs, where the
  long parameter list is followed by the brace on a separate line;
- one space after `if`, `for`, `while`, `switch`, none after the function
  name in calls;
- the pointer star is attached to the name (`char *s`), also in return
  types where readable (`Logger* new_logger(...)` is accepted);
- comments: `/* ... */` for documentation, `//` for short notes inside
  functions;
- function-like macros which are statements are wrapped in
  `do { ... } while (0)`, so that they can be used anywhere a statement
  can (e.g. in an `if` without braces); their arguments are enclosed in
  parentheses in the expansion;
- functions which never return (e.g. that report an error and exit) are
  marked `__attribute__ ((noreturn))`, so that the compiler and the static
  analyzers know it;
- each `.c` file includes its own header, so that the declarations are
  checked against the definitions; functions without parameters are
  declared with `(void)`.

### Naming

- functions and variables: `lower_snake_case`;
- macros and constants: `UPPER_SNAKE_CASE`, prefixed by the module name
  (e.g. `READ_BATCH_INIT_SIZE`);
- types: `module_name_t` (e.g. `samples_info_t`); enumeration values in
  upper case with a common prefix (e.g. `LOG_ERROR`, `LOG_WARN`);
- functions of a module are prefixed by the module name:
  `samples_info_new()`, `samples_info_match()`, `samples_info_delete()`;
- internal (non-exported) functions are `static` and do not need the prefix;
  those which must be in a header (used by `static inline` functions) end
  with `_` (e.g. `short_edit_rows_()`), and are not part of the API;
- the object on which a function operates is the first parameter and is
  called `self`;
- functions which are wrapped by a macro adding the source location are
  prefixed by `_` (see below); names starting with `__` are reserved and
  not used.

### Includes

Standard library headers first, then the project headers, each group in
any order; include what is used, do not rely on indirect includes.

## Modules

A module encapsulates one data structure or one task. The structure is
opaque: the header only declares the type, and its fields are defined in
the `.c` file.

Header:

```c
typedef struct index_table_t index_table_t;

/*
   Create a new table.

   Input parameters:
   - <n_indices> number of indices for each sample
   - <logger> logger for the verbose messages (may be NULL)
*/
index_table_t* index_table_new(size_t n_indices, Logger *logger);

/* Add the indices of a sample; the strings are copied */
void index_table_add(index_table_t *self, const char *sample_id,
    char *const *indices);

/* Delete a table (NULL is accepted) */
void index_table_delete(index_table_t *self);

/* Unit tests */
void index_table_test(void);
```

Implementation:

```c
#define INDEX_TABLE_INIT_SIZE 64
#define INDEX_TABLE_INC_SIZE  64

struct index_table_t {
  index_table_item_t *items;
  size_t len;
  size_t cap;
};

index_table_t* index_table_new(size_t n_indices, Logger *logger) {
  index_table_t *self;
  XMALLOC(self, sizeof(index_table_t));
  self->len = 0;
  self->cap = INDEX_TABLE_INIT_SIZE;
  XMALLOC(self->items, sizeof(index_table_item_t) * self->cap);
  ...
  return self;
}

void index_table_add(index_table_t *self, const char *sample_id,
    char *const *indices) {
  assert(self != NULL);
  assert(sample_id != NULL);
  ...
}
```

Rules:

- each type has a constructor `_new()` and a destructor `_delete()`; the
  destructor accepts NULL and frees everything the object owns;
- every exported function asserts its preconditions (non-NULL pointers,
  valid ranges);
- every function is documented in the header, before its declaration:
  what it does, the meaning of the parameters (`<name>` notation), the
  ownership of the returned values (who must free them) and the behaviour
  on errors;
- no global variables, except constants; state is kept in objects, so
  that modules can be used from several threads, each with its own
  objects;
- dynamic arrays are represented by a pointer plus `len` and `cap` fields
  and grown with the `XCHECKSPACE*` macros;
- exception to the opacity: plain data structures, which are only read
  after their construction (e.g. a configuration), can be defined in the
  header, so that their fields are accessed directly.

## Memory allocation

Memory is allocated only through the macros of `alloc.h`, which exit the
program with an error message if the allocation fails (so that the return
value never needs to be checked):

- `XMALLOC(PTR, SIZE)`, `XCALLOC(PTR, NITEMS, ITEMSIZE)`,
  `XREALLOC(PTR, SIZE)`;
- `XFREE(PTR)`: frees the memory and sets the pointer to NULL;
- `XCHECKSPACE(PTR, CAPPTR, REQSIZE)`: grows a buffer to at least the
  required size; `XCHECKSPACEINC` grows it to the required size plus a
  constant increment, `XCHECKSPACEINCPERC` plus a percentage of the current
  size (to obtain amortized constant time appends).

The sizes passed to these macros are in bytes.

In the inner loops, executed millions of times (e.g. the processing of
each record of a large file), no memory is allocated per iteration:
buffers are allocated once (per thread) and reused, growing them when
needed. The memory usage must not depend on the size of the input.

## Logging and errors

All messages are output through a `Logger` (`logger.h`), created by the
program from the command line options and passed to the modules that need
it. Messages go to the standard error (or to a log file); the standard
output is reserved for the results.

Levels, from the most to the least important:

| level | meaning |
|---|---|
| `LOG_ERROR` | something went wrong, the program is interrupted |
| `LOG_WARN` | something went wrong, but the program continues |
| `LOG_SUCCESS` | an operation was completed successfully |
| `LOG_PHASE` | start and end of a phase of the program |
| `LOG_INFO` | additional information (verbose mode) |
| `LOG_PARAM` | the values of the parameters (verbose mode) |
| `LOG_MEM` | memory usage (verbose mode) |
| `LOG_DEBUG` | debugging information, with low impact on speed |
| `LOG_TRACE` | debugging information, with high impact on speed |

When the output is a terminal, the levels are shown in different colors;
each message is preceded by the time and the level.

Usage:

- `lprintf(logger, LOG_INFO, "Number of samples: %zu\n", n)`: output a
  message (`lnprintf`: without newline); the macros add the source file and
  line, shown at the debug levels;
- `LCHECK(condition, logger, "message", ...)`: if the condition is false,
  log an error and exit; used for checking the input and the results of
  system calls;
- `fatal_printf("message", ...)`: exit with an error where no logger is
  available;
- `showmem(logger)`: log the current memory usage;
- `TRACE`/`NTRACE` and, compiled only with `-DSHOW_TRACES`, `DEBUG_TRACE`:
  temporary debugging output.

The modules never write messages directly to the standard error (with
`fprintf(stderr, ...)`) and do not call `exit()` without logging the
error; an error which is not reported through `LCHECK` is logged with
`lprintf(logger, LOG_ERROR, ...)`, followed by `exit(EXIT_FAILURE)`.
The level is shown by the logger, so the messages do not start with
"Error". The functions which open files with a library other than the C
standard library (e.g. zlib) report the errors in the same way.

The logger argument can always be NULL, in which case only errors are
output. Error messages describe the problem in terms of the user's input
(file name, line, value), not of the implementation.

Errors in the input lead to an exit with `EXIT_FAILURE`; the functions thus
do not return error codes, unless the caller is expected to handle the
error (in which case this is documented in the header).

## File operations

Files are opened, read, written and closed through the wrappers of `file.h`,
which have the name of the standard function with a prepended `x`
(`xfopen`, `xfclose`, `xfread`, `xfwrite`, `xgetline`, `xfunlink`,
`xfrename`, ...). Differences to the standard functions:

- on error, the program is interrupted with a descriptive message
  (including the file name);
- the file name suffix is passed separately (or NULL), so that the caller
  does not need to concatenate strings; the suffixes are constants defined
  in the `filesfx.h` of each project, starting with a dot (e.g.
  `FILESFX_SORTED`);
- a description of the file (used in the messages) and a logger are
  passed;
- the implementation function is called `_xfopen` and takes the
  source file and line, added by the `xfopen` macro.

The same pattern (a `_` function and a macro adding `__FILE__` and
`__LINE__`) is used for any other function that reports errors on behalf
of its caller, e.g. `_lprintf` and `_LCHECK` of the logger, behind the
`lprintf` and `LCHECK` macros.

## HTML output

HTML (and inline SVG) is not written as markup in C string literals, but:

- static markup (page skeletons, fixed blocks of text) as templates:
  HTML files, indented normally, with named placeholders (`{{name}}`,
  escaped; `{{{name}}}`, trusted markup; `{{#name}}...{{/name}}`,
  optional sections), embedded at compile time (see "Building") and
  written by `html_template()` (`html-template.h`);
- dynamic markup (tables, plots) by the element writer `html_writer_t`
  (`html-writer.h`): one call to open an element, one for each attribute
  (with a printf format checked by the compiler), one to close it; the
  calls are indented in the code as the elements in the page:

```c
html_open(w, "tr");
  html_inline(w);
  html_tag_text(w, "td", name);
  html_open(w, "td");
    html_textf(w, "%s", html_num(count).s);
  html_close(w);
html_close(w);
```

Text and attribute values are escaped by the writer (or the template), so
that no other code escapes HTML; trusted markup is passed explicitly with
`html_raw()` or `{{{name}}}`. The output is indented, one element per
line (elements containing text on one line), so that a report can be read
and compared with a diff.

## Command line programs

Each program is a `src/<program>.c` file with the following parts:

1. the includes and the default values of the options
   (`#define DEFAULT_MAX_ERRORS 0.25`);
2. the help message, as a `PURPOSE` macro, marked by a `// %-PURPOSE-%`
   comment, with a `%s` for the program name:

   ```c
   // %-PURPOSE-%
   #define PURPOSE "Extract the inserts from the reads\n"\
     "\n"\
     "Usage:\n"\
     "   %s [options] <reads.fastq[.gz]>\n"\
     "\n"\
     "Arguments:\n"\
     "   <reads.fastq[.gz]>  input reads\n"\
     "\n"\
     "Options:\n"\
     "   -t <n>  number of threads (default: " XSTR(DEFAULT_THREADS) ")\n"\
     "\n"\
     COMMON_OPTS_HELP
   ```

3. a `parseargs()` function, which sets the output parameters to their
   defaults, parses the options with `getopt()`, checks the number of
   arguments, creates the logger and logs the parameter values:

   ```c
   void parseargs(int argc, char *argv[],
       char **reads_filename,
       size_t *threads,
       Logger **logger)
   {
     PARSE_START;
     int c = 0;
     *threads = DEFAULT_THREADS;
     opterr = 0;
     while ((c = getopt(argc, argv, COMMON_OPTS_STR "t:")) != -1) {
       switch (c) {
         case 't':
           PARSE_INTEGER_OPT_MIN('t', threads, 1);
           break;
         COMMON_OPTS_PARSE_CASES;
       }
     }
     CHECK_NARGS(1);
     *reads_filename = argv[optind];
     INIT_LOGGER(*logger);
     PARSE_LOG("Reads file: %s", *reads_filename);
     PARSE_LOG("Threads: %zu", *threads);
   }
   ```

4. a short `main()`, which calls `parseargs()`, runs the work by calling
   the modules, frees all memory and returns `EXIT_SUCCESS`.

The macros of `argparse-def.h` give all programs the same interface:

- common options: `-q` (quiet), `-v` (verbose), `-D` (debug),
  `-L <file>` (log file), `-h` (help) and `-V` (version);
- `OPT_ERR(C, ...)`, `ARG_ERR(N, ...)`, `PARSE_ERR(...)`: report an error
  in an option or argument, followed by the usage message;
- `PARSE_INTEGER_OPT_MIN`, `PARSE_INTEGER_OPT_MIN_MAX`,
  `PARSE_INTEGER_ARG_MIN`: parse and check integer values;
- `CHECK_NARGS(N)`, `CHECK_MIN_NARGS(N)`: check the number of positional
  arguments;
- `VERSION`, `STR`/`XSTR` (stringify), `BOOLSTR`.

Options are single letters, parsed with the standard `getopt()`.

## Tests

### Unit tests

Modules with non-trivial logic export a `<module>_test()` function, which
uses the macros of `test.h`:

- `test_assert(condition)`: exit with an error message if false;
- `test_eq(found, expected, type, format)` and the shortcuts
  `test_eq_size`, `test_eq_int`, `test_eq_long`, `test_eq_char`,
  `test_eq_str`: compare two values and output `[OK]` or the differing
  values;
- `test_eq_mem(found, len, expected)`: the same for `len` characters which
  are not 0-terminated (e.g. the content of a `strbuf_t`).

The unit tests are run by a test program (e.g. `bin/cbaia-tests`), with one
option for each test group (listed in a table), so that they can be run
separately; without options, all tests are run. The options must not
collide with the common options (e.g. `-q`). The test program is run
under valgrind, failing on memory errors and leaks.

### Program tests

The programs are tested on small input files in `testdata/`, comparing their
output with the expected results (files with suffix `.expected`), through the
`run_cmp` function of `scripts/test-cmd-lib.sh`:

- the program is run under valgrind (`--error-exitcode=1 --leak-check=full
  --errors-for-leak-kinds=definite`), so that memory errors and leaks make
  the test fail;
- the exit code is checked;
- up to three output files are compared with `diff` (optionally after
  sorting them); on differences, the output is copied to
  `<expected>.found`, for inspection.

Each test is a Make target (`test_<group>_<case>`); the targets of a group
are generated from lists of parameters with `foreach`, and the failing
targets are listed at the end of `make test`, so that they can be repeated
individually.

When a program reimplements an existing tool (e.g. a prototype in
another language), the expected results are produced by that tool: the C
programs must produce byte-identical outputs on the same inputs and
parameters.

### Other checks

Before a release, and after changes to the parsers or to the memory
handling:

- compile and run the tests with `-fsanitize=address,undefined`;
- compile with `-Wall -Wextra` without warnings (e.g. in the continuous
  integration with `-Werror`);
- run a static analyzer, e.g. `clang-tidy
  --checks='clang-analyzer-*,bugprone-*'`, and review its findings;
- fuzz the input parsers (`make fuzz`, then `tests/fuzz/run.sh <harness>
  <seconds> [<jobs>]`; harnesses: `reads`, `bam`, `whitelist`); a new
  parser comes with a harness in `tests/fuzz`.

## Performance

- Measure before optimizing: use the `DEBUGFLAGS_PROF` build, or `perf`, on
  realistic inputs; record the benchmark results with the commit used.
- Keep the inner loops free of allocations, I/O and logging (use
  `LOG_TRACE` only in debug builds).
- Prefer simple, cache-friendly data structures (arrays with `len`/`cap`)
  to linked structures.
- Multithreading: the input is read in batches (under a lock), each
  thread processes its batch with its own objects (no shared mutable
  state), and the results are written in the input order (each thread
  waits for its turn), so that the output does not depend on the number of
  threads.

## Commits

- Commits are atomic: one logical change for each commit (e.g. a new
  module with its unit tests; a bug fix with its test), and each commit
  compiles and passes the tests.
- The commit message has a short summary line (imperative mood), a blank
  line and, if needed, a description of the reason of the change.
- Everything is developed on the `main` branch.
