/* (c) 2026 Giorgio Gonnella */

/*
   The built-in reference set: complete segments (with the UTRs) of
   several influenza A subtypes and of influenza B, so that a run can
   be classified against references without the user providing any. The
   set is not shipped with flucut; it is downloaded from NCBI on first
   use and cached, as it is the same for every run and may be updated
   independently of the program.

   builtin_references_path() returns the path to the cached FASTA of
   the types <set> (--virus, see virus.h), downloading it if it is not
   there yet; NULL if it cannot be obtained (no network, no curl): the
   caller then runs without references. Each set of types has its own
   cache file, so that a run of one type does not download the
   references of the other.

   The environment variable FLUCUT_BUILTIN_REFERENCES overrides the set
   with a local file (no download), for offline use and the tests, and
   the key "references" of the configuration file names a set to use
   instead of it (see user-config.h), for instance one built from the
   seed sequences of IRMA.
*/

#ifndef BUILTIN_REFERENCES_H
#define BUILTIN_REFERENCES_H

#include "logger.h"
#include "virus.h"

const char* builtin_references_path(virus_set_t set, Logger *logger);

#endif
