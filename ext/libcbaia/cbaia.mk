#
# Makefile fragment for the projects using a copy of libcbaia (vendoring).
#
# Usage, in the Makefile of the project:
#
#   CBAIA_DIR = ext/libcbaia
#   include $(CBAIA_DIR)/cbaia.mk
#
#   CFLAGS += $(CBAIA_INCLUDES)
#   LDLIBS += $(CBAIA_LDLIBS)
#   bin/prog: src/prog.c $(CBAIA_LIBS) ...
#   	$(CC) ... $(CBAIA_LIB) $(LDLIBS)
#
# The library is built (by the Makefile of libcbaia, with the DEBUGFLAGS
# of the project, if defined) when its sources change.
#
# zlib: by default, the copy included in libcbaia is used; with
# SYSTEM_LIBS=1 (e.g. make SYSTEM_LIBS=1), the system library.
#

# the components of the library (subdirectories of src), also used by
# its Makefile
CBAIA_COMPONENTS = core ds io seq
CBAIA_INCLUDES = $(addprefix -I$(CBAIA_DIR)/src/,$(CBAIA_COMPONENTS)) \
  -I$(CBAIA_DIR)/ext/khash
CBAIA_LIB = $(CBAIA_DIR)/lib/libcbaia.a
ifeq ($(SYSTEM_LIBS),1)
CBAIA_LDLIBS = -lz -lm
CBAIA_ZLIB =
else
CBAIA_INCLUDES += -I$(CBAIA_DIR)/ext/zlib
CBAIA_ZLIB = $(CBAIA_DIR)/lib/libcbaia-zlib.a
CBAIA_LDLIBS = $(CBAIA_ZLIB) -lm
endif
# all libraries built by libcbaia (prerequisites of the programs)
CBAIA_LIBS = $(CBAIA_LIB) $(CBAIA_ZLIB)

CBAIA_MAKE_LIB = $(MAKE) -C $(CBAIA_DIR) lib SYSTEM_LIBS=$(SYSTEM_LIBS) \
  $(if $(DEBUGFLAGS),DEBUGFLAGS='$(DEBUGFLAGS)')

# the Makefile of the library is always run: it rebuilds the library if
# its sources or the compilation flags (DEBUGFLAGS) changed; otherwise the
# library is not modified, and the programs are not relinked
$(CBAIA_LIB): cbaia-force
	$(CBAIA_MAKE_LIB)

.PHONY: cbaia-force
cbaia-force:

# the zlib library is built by the same run of the Makefile of the library
# (with make -j, two runs at the same time would conflict); its recipe is
# also needed when only it is missing (e.g. after a build with
# SYSTEM_LIBS=1)
ifneq ($(CBAIA_ZLIB),)
$(CBAIA_ZLIB): | $(CBAIA_LIB)
	$(CBAIA_MAKE_LIB)
endif

# script writing a data file (e.g. an HTML template) as a C string
# literal, see dev/embed.sh; e.g. for the files src/*.html of the project:
#   obj/%-html.inc: src/%.html $(CBAIA_EMBED)
#   	$(CBAIA_EMBED) $< > $@
# with .DELETE_ON_ERROR: in the Makefile of the project, so that a failure
# does not leave a truncated .inc, considered up to date by the next run;
# it is not set here, as it applies to all the rules of the Makefile
# including this fragment (the library itself is built by its own
# Makefile, which sets it)
CBAIA_EMBED = $(CBAIA_DIR)/dev/embed.sh

.PHONY: cbaia-clean
cbaia-clean:
	$(MAKE) -C $(CBAIA_DIR) cleanup
