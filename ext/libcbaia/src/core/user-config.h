/* (c) 2026 Giorgio Gonnella */

/*
   The configuration file of flucut: $XDG_CONFIG_HOME/flucut/config, or
   ~/.config/flucut/config if XDG_CACHE_HOME is not set. It holds lines
   of "key = value"; a # starts a comment, and keys the program does not
   know are ignored, so that the tools around flucut can keep their own
   settings in the same file (scripts/irma_references.py writes where
   IRMA and its seed sequences are, to build the reference set again
   when IRMA updates them).

   The program itself reads one key, "references": the reference set of
   a run without -r, used instead of the built-in set (see
   builtin-references.h). scripts/irma_references.py writes it when it
   builds a reference set from the seeds of IRMA.
*/

#ifndef USER_CONFIG_H
#define USER_CONFIG_H

#include "logger.h"

/*
   The value of <key> in $XDG_CONFIG_HOME/<program>/config or
   ~/.config/<program>/config, or NULL if the file or key is not there.
   The result is in a static buffer, valid until the next call.
*/
const char* user_config_value(const char *program, const char *key,
    Logger *logger);

/* Unit tests */
void user_config_test(void);

#endif
