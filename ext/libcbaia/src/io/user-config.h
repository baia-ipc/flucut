/* (c) 2026 Giorgio Gonnella */

/*
   A program-specific configuration file at
   $XDG_CONFIG_HOME/<program>/config, or ~/.config/<program>/config if
   XDG_CONFIG_HOME is not set. It holds lines of "key = value"; a #
   starts a comment, and unknown keys are ignored so several tools can
   share a configuration file.
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
