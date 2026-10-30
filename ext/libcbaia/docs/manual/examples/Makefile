#
# Example programs of the libcbaia manual, built with the copy of the
# library of this repository (as a project vendoring it would do).
#
#   make         build the examples
#   make check   build and run them on the data in data/, comparing the
#                output with the expected results
#

PROGS = readcount bcdemux logdemo htmlreport colors

# the default goal (before the rules of cbaia.mk)
all: $(PROGS)

CBAIA_DIR = ../../..
include $(CBAIA_DIR)/cbaia.mk

CFLAGS += -Wall -Wextra -std=gnu11 -O2 $(CBAIA_INCLUDES)
LDLIBS += $(CBAIA_LDLIBS)

%: %.c $(CBAIA_LIBS)
	$(CC) $(CFLAGS) -o $@ $< $(CBAIA_LIB) $(LDLIBS)

VALGRIND = $(if $(shell command -v valgrind),valgrind --error-exitcode=1 -q \
           --leak-check=full --errors-for-leak-kinds=definite)

check: $(PROGS)
	rm -rf out && mkdir out
	$(VALGRIND) ./readcount -q data > out/readcount.tsv
	diff out/readcount.tsv expected/readcount.tsv
	$(VALGRIND) ./bcdemux -q -w data/barcodes.txt -d 1 -o out \
	  data/reads.fastq
	for f in expected/reads.*.fastq; do \
	  zcat out/`basename $$f`.gz | diff - $$f || exit 1; done
	$(VALGRIND) ./logdemo -n 1000 2> out/logdemo.err
	sed 's/^# [0-9-]* [0-9:]* |/#/' out/logdemo.err | diff - expected/logdemo.log
	! ./logdemo -q -n 10 -e 2> /dev/null
	$(VALGRIND) ./htmlreport -q out/report.html
	grep -v '>Updated: ' out/report.html | diff - expected/report.html
	$(VALGRIND) ./colors > out/colors.txt
	diff out/colors.txt expected/colors.txt
	rm -rf out
	@echo "[OK] examples"

clean:
	rm -rf $(PROGS) out

.PHONY: all check clean
