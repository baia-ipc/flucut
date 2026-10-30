/* (c) 2026 Giorgio Gonnella */

/*
   The banner of flucut: the logo as colored text (Unicode Braille
   characters and 24-bit ANSI colors, 77 columns, 13 lines), shown at the
   start of a run in a terminal. Its source is docs/logo/flucut-banner.ans
   (see docs/logo/README.md), embedded at compile time.
*/

#ifndef BANNER_H
#define BANNER_H

#include <stdio.h>

/* Write the banner to <file> */
void banner_print(FILE *file);

#endif
