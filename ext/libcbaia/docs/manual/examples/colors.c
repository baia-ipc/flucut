/* (c) 2026 Giorgio Gonnella */

/*
   Example program of the libcbaia manual: the foreground colors of
   colors.h, each name written in its color.
*/

#include <stdio.h>
#include "colors.h"

int main(void) {
  static const struct { const char *code, *name; } colors[] = {
    {COLOR_BLACK, "COLOR_BLACK"}, {COLOR_RED, "COLOR_RED"},
    {COLOR_GREEN, "COLOR_GREEN"}, {COLOR_YELLOW, "COLOR_YELLOW"},
    {COLOR_BLUE, "COLOR_BLUE"}, {COLOR_MAGENTA, "COLOR_MAGENTA"},
    {COLOR_CYAN, "COLOR_CYAN"}, {COLOR_WHITE, "COLOR_WHITE"},
    {COLOR_BRIGHT_BLACK, "COLOR_BRIGHT_BLACK"},
    {COLOR_BRIGHT_RED, "COLOR_BRIGHT_RED"},
    {COLOR_BRIGHT_GREEN, "COLOR_BRIGHT_GREEN"},
    {COLOR_BRIGHT_YELLOW, "COLOR_BRIGHT_YELLOW"},
    {COLOR_BRIGHT_BLUE, "COLOR_BRIGHT_BLUE"},
    {COLOR_BRIGHT_MAGENTA, "COLOR_BRIGHT_MAGENTA"},
    {COLOR_BRIGHT_CYAN, "COLOR_BRIGHT_CYAN"},
    {COLOR_BRIGHT_WHITE, "COLOR_BRIGHT_WHITE"},
  };
  for (size_t i = 0; i < sizeof(colors) / sizeof(colors[0]); i++) {
    START_COLOR(stdout, colors[i].code);
    printf("%s", colors[i].name);
    STOP_COLOR(stdout);
    printf("\n");
  }
  return 0;
}
