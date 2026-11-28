/* (c) 2026 Giorgio Gonnella */

#include <stdio.h>
#include "banner.h"

static const char *const banner =
#include "flucut-banner-ans.inc"
  ;

void banner_print(FILE *file) {
  fputs(banner, file);
}
