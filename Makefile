BINDIR=bin
OBJDIR=obj

SRCDIR=src
SRCS=$(wildcard $(SRCDIR)/*.c)

PROGS = $(BINDIR)/flucut $(BINDIR)/fc-tests

.PHONY: default objs install uninstall clean cleanup cleandep cbaia-clean \
        test check_style test_units test_programs test_options test_live \
        test_scripts \
        manual

# (a failed rule, e.g. the embedding of a template, leaves no partial file)
.DELETE_ON_ERROR:

default: $(PROGS)

EXTDIR=ext

CBAIA_DIR = $(EXTDIR)/libcbaia
include $(CBAIA_DIR)/cbaia.mk

# libyaml: vendored (default) or the system library (make SYSTEM_LIBS=1);
# the same holds for zlib, included in libcbaia (see cbaia.mk)
YAMLDIR = $(EXTDIR)/libyaml
ifeq ($(SYSTEM_LIBS),1)
YAML_INCLUDES =
YAML_OBJS =
YAML_LDLIBS = -lyaml
else
YAML_INCLUDES = -I$(YAMLDIR)
YAML_OBJS = $(patsubst $(YAMLDIR)/%.c,$(OBJDIR)/yaml/%.o,\
              $(wildcard $(YAMLDIR)/*.c))
YAML_LDLIBS =
endif

# libdeflate (compression of the output): always the vendored copy, also
# with SYSTEM_LIBS=1 (it is not installed on most systems); all the
# sources of the library, including those of the other architectures
# (their content is compiled only on them)
DEFLATEDIR = $(EXTDIR)/libdeflate
DEFLATE_OBJS = $(patsubst $(DEFLATEDIR)/%.c,$(OBJDIR)/libdeflate/%.o,\
                 $(wildcard $(DEFLATEDIR)/lib/*.c $(DEFLATEDIR)/lib/*/*.c))

CFLAGS += -Wall -Wextra -std=gnu11 -I$(EXTDIR) -I$(EXTDIR)/edlib \
          -I$(DEFLATEDIR) $(CBAIA_INCLUDES) $(YAML_INCLUDES)
CXXFLAGS += -std=c++11 -I$(EXTDIR)/edlib
LDLIBS += $(YAML_LDLIBS) $(CBAIA_LDLIBS) -lpthread

# external libraries, linked into all programs
EXTOBJS = $(OBJDIR)/edlib.o $(YAML_OBJS) $(DEFLATE_OBJS) $(CBAIA_LIBS)

DEBUGFLAGS_DEBUG = -gdwarf-3 -g3 -O0
DEBUGFLAGS_NDEBUG = -O3 -DNDEBUG -flto
DEBUGFLAGS_PROF = -no-pie -pg $(DEBUGFLAGS_DEBUG)
DEBUGFLAGS_SANITIZE = -fsanitize=address,undefined $(DEBUGFLAGS_DEBUG)
# select the DEBUGFLAGS settings to use:
DEBUGFLAGS=$(DEBUGFLAGS_NDEBUG)

# header dependency tracking recipy modified from:
# http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/

DEPDIR := .d
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td

COMPILE.c = $(CC) $(DEPFLAGS) $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(DEBUGFLAGS) -c
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

$(OBJDIR)/%.o : $(SRCDIR)/%.c
$(OBJDIR)/%.o : $(SRCDIR)/%.c $(DEPDIR)/%.d | $(OBJDIR) $(DEPDIR)
	$(COMPILE.c) $(OUTPUT_OPTION) $<
	$(POSTCOMPILE)

# the HTML templates of the report (src/*.html), embedded in the program:
# src/<name>.html gives obj/<name>-html.inc (see $(CBAIA_EMBED)); the
# objects including them depend on them through the header dependency
# tracking, after the first build
TEMPLATE_INCS = $(patsubst $(SRCDIR)/%.html,$(OBJDIR)/%-html.inc,\
                  $(wildcard $(SRCDIR)/*.html))
CFLAGS += -I$(OBJDIR)
$(OBJDIR)/%-html.inc: $(SRCDIR)/%.html $(CBAIA_EMBED) | $(OBJDIR)
	$(CBAIA_EMBED) $< > $@
$(OBJDIR)/html-report.o: | $(TEMPLATE_INCS)

# the banner shown in a terminal (docs/logo/flucut-banner.ans), embedded
# in the same way
BANNER_INC = $(OBJDIR)/flucut-banner-ans.inc
$(BANNER_INC): docs/logo/flucut-banner.ans $(CBAIA_EMBED) | $(OBJDIR)
	$(CBAIA_EMBED) $< > $@
$(OBJDIR)/banner.o: | $(BANNER_INC)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(notdir $(basename $(SRCS)))))

# end of header dependency tracking

# the output directories, created when needed (order-only prerequisites:
# their modification times are ignored)
$(BINDIR) $(OBJDIR) $(OBJDIR)/yaml $(DEPDIR):
	mkdir -p $@

$(OBJDIR)/edlib.o: $(EXTDIR)/edlib/edlib.cpp $(EXTDIR)/edlib/edlib.h \
    | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(DEBUGFLAGS) -c -o $@ $<

YAML_CFLAGS = -O2 -DYAML_VERSION_MAJOR=0 -DYAML_VERSION_MINOR=2 \
              -DYAML_VERSION_PATCH=5 -DYAML_VERSION_STRING='"0.2.5"'

$(OBJDIR)/yaml/%.o: $(YAMLDIR)/%.c $(YAMLDIR)/yaml.h $(YAMLDIR)/yaml_private.h \
    | $(OBJDIR)/yaml
	$(CC) $(YAML_CFLAGS) -I$(YAMLDIR) -c -o $@ $<

# (the flags of the build system of libdeflate, for a release build; the
# objects are in subdirectories, as the sources: lib/x86, lib/arm)
DEFLATE_CFLAGS = -O2 -I$(DEFLATEDIR)
DEFLATE_HEADERS = $(DEFLATEDIR)/libdeflate.h $(DEFLATEDIR)/common_defs.h \
                  $(wildcard $(DEFLATEDIR)/lib/*.h $(DEFLATEDIR)/lib/*/*.h)

$(OBJDIR)/libdeflate/%.o: $(DEFLATEDIR)/%.c $(DEFLATE_HEADERS)
	@mkdir -p $(@D)
	$(CC) $(DEFLATE_CFLAGS) -c -o $@ $<

# objects to link into each program, computed from the included headers
# (see dev/objdeps.sh of libcbaia)

# (all objects are built first, by objs, not by a sub-make for each
# program: with make -j, the sub-makes would build the same objects at the
# same time; with an empty list, a sub-make would build the default goal)
$(DEPDIR)/%.o_list: $(SRCDIR)/%.c $(SRCS) $(wildcard $(SRCDIR)/*.h) \
    | objs $(DEPDIR)
	$(CBAIA_DIR)/dev/objdeps.sh $(SRCDIR) $(OBJDIR) $< > $@

# linked by the C++ compiler, for the C++ standard library used by edlib;
# the object of the program (with main()) is built by the rule of all
# objects (thus with the header dependency tracking), and is not in the
# list of objects (there is no header <program>.h); relinked also when
# an object changes (e.g. after a change of a header of the vendored
# libcbaia, which does not change the list of objects)
$(BINDIR)/%: $(OBJDIR)/%.o $(DEPDIR)/%.o_list $(EXTOBJS) \
    $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SRCS)) | $(BINDIR)
	$(CXX) $(DEBUGFLAGS) $(LDFLAGS) -o $@ $(OBJDIR)/$*.o \
	  `cat $(DEPDIR)/$*.o_list` $(EXTOBJS) $(LDLIBS)

objs: $(EXTOBJS) $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SRCS))

PREFIX ?= /usr/local

install: $(BINDIR)/flucut
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	cp $(BINDIR)/flucut $(DESTDIR)$(PREFIX)/bin/

uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/flucut

cleanup: clean cbaia-clean
	rm -rf bin

clean: cleandep
	rm -rf obj

cleandep:
	rm -rf .d

# the tests are run under valgrind, if available (e.g. not on recent macOS)
VALGRIND_EXITCODE = 1
VALGRIND = $(if $(shell command -v valgrind),valgrind -q \
           --error-exitcode=$(VALGRIND_EXITCODE) \
           --leak-check=full --errors-for-leak-kinds=definite)

test: check_style $(PROGS) test_units test_programs test_options \
      test_live test_scripts

# formatting rules of the developer guide (line length, whitespace, file
# headers, include guards, indentation), for the C sources, the HTML
# templates and the Python and shell scripts
check_style:
	@tests/check-style.sh

test_units: $(BINDIR)/fc-tests
	$(VALGRIND) $(BINDIR)/fc-tests -q

# program tests: the cases of testdata/cases.tsv, run under valgrind; the
# outputs are compared to the expected outputs (testdata/expected); on
# differences, the output is written to
# testdata/expected/<case>.txt.found
CASES := $(shell grep -v '^\#' testdata/cases.tsv | cut -f1)

TEST_PROGRAMS := $(addprefix test_program_,$(CASES))
.PHONY: $(TEST_PROGRAMS)

test_programs: $(TEST_PROGRAMS)

# validation of the options: invalid values are rejected (exit status 1,
# thus another one for the errors found by valgrind)
test_options: VALGRIND_EXITCODE = 99
test_options: $(BINDIR)/flucut
	@tests/test-options.sh $(BINDIR)/flucut "$(VALGRIND)"

# live mode: files written during the run, whose outputs must be
# identical to those of a normal run on the final input; resume,
# interruption, end of the run, -e auto (tests/test-live.sh)
test_live: $(BINDIR)/flucut
	@tests/test-live.sh $(BINDIR)/flucut "$(VALGRIND)"

# the scripts, offline (the download of the references from NCBI:
# python3 -I tests/test-fetch-references.py --live)
test_scripts:
	@python3 -I tests/test-fetch-references.py

# (a static pattern rule: the implicit rules are not searched for the
# phony targets)
$(TEST_PROGRAMS): test_program_%: $(BINDIR)/flucut
	@tests/run-case.sh $(BINDIR)/flucut $* "$(VALGRIND)" \
	  > testdata/expected/$*.txt.found
	@if cmp -s testdata/expected/$*.txt.found testdata/expected/$*.txt; \
	  then rm testdata/expected/$*.txt.found; echo "[OK] program test $*"; \
	  else echo "[FAIL] program test $*: see" \
	    "testdata/expected/$*.txt.found"; exit 1; fi

# user manual (PDF), from docs/manual (requires Sphinx and LaTeX)
SPHINXBUILD ?= sphinx-build
manual:
	$(SPHINXBUILD) -M latexpdf docs/manual docs/manual/_build -q
	mkdir -p manual
	cp docs/manual/_build/latex/flucut.pdf manual/flucut-manual.pdf
