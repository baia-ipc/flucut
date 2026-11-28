/* (c) 2026 Giorgio Gonnella */

/*
   The adapters at the ends of a read without a primer (see "Adapters" in
   amplicon.h): a step of the analysis of a read (amplicon-analyzer.h).

   The read is analyzed; at an end without a primer the adapters are
   searched (adapters.h) and, if one is found, the part between the
   adapters is analyzed again, with the primer of an end next to its
   adapter within AMPLICON_ADAPTER_WINDOW bases. The second analysis is
   taken if it finds at least as many primers as the first.
*/

#ifndef AMPLICON_ADAPTERS_H
#define AMPLICON_ADAPTERS_H

#include <stddef.h>
#include "amplicon.h"

/* The analysis of the read <seq> (length <len>), with the adapters, into
   <out> */
void amplicon_analyze_with_adapters(amplicon_analyzer_t *self,
    const char *seq, size_t len, amplicon_t *out);

#endif
