/* (c) 2026 Giorgio Gonnella */

/*
   The built-in reference set: complete segments (with the UTRs) of several
   influenza A subtypes, so that a run can be classified against references
   without the user providing any. The set is not shipped with flucut; it
   is downloaded from NCBI on first use and cached, as it is the same for
   every run and may be updated independently of the program.

   builtin_references_path() returns the path to the cached FASTA,
   downloading it if it is not there yet; NULL if it cannot be obtained
   (no network, no curl): the caller then runs without references.

   The environment variable FLUCUT_BUILTIN_REFERENCES overrides the set
   with a local file (no download), for offline use and the tests.
*/

#ifndef BUILTIN_REFERENCES_H
#define BUILTIN_REFERENCES_H

#include "logger.h"

const char* builtin_references_path(Logger *logger);

#endif
