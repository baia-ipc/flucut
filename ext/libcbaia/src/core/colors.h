/* ANSI colors of the text written to a terminal */

#ifndef COLORS_H
#define COLORS_H

#include <stdio.h>

/* write the text after it in COLOR (one of the COLOR_ codes below) */
#define START_COLOR(STREAM, COLOR) \
  fprintf(STREAM, "\x1b[%sm", COLOR)

/* back to the default color */
#define STOP_COLOR(STREAM) \
  fprintf(STREAM, "\x1b[0m")

#define COLOR_BLACK          "30"
#define COLOR_RED            "31"
#define COLOR_GREEN          "32"
#define COLOR_YELLOW         "33"
#define COLOR_BLUE           "34"
#define COLOR_MAGENTA        "35"
#define COLOR_CYAN           "36"
#define COLOR_WHITE          "37"

#define COLOR_BRIGHT_BLACK   "90"
#define COLOR_BRIGHT_RED     "91"
#define COLOR_BRIGHT_GREEN   "92"
#define COLOR_BRIGHT_YELLOW  "93"
#define COLOR_BRIGHT_BLUE    "94"
#define COLOR_BRIGHT_MAGENTA "95"
#define COLOR_BRIGHT_CYAN    "96"
#define COLOR_BRIGHT_WHITE   "97"

#endif
