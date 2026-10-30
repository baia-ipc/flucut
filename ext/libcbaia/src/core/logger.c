/* (c) 2021 Giorgio Gonnella */

#include <errno.h>
#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <sys/resource.h>
#include "alloc.h"
#include "colors.h"
#include "logger.h"

struct Logger {
  FILE *logfile;
  LoggingLevel level;
  bool colors;
};

Logger* new_logger(FILE *logfile, LoggingLevel level) {
  Logger *result;
  XMALLOC(result, sizeof(*result));
  result->logfile = logfile;
  result->level = level;
  result->colors = isatty(fileno(logfile));
  return result;
}

void delete_logger(Logger* logger) {
  if (logger != NULL) {
    if (logger->logfile != stderr && logger->logfile != stdout) {
      if (fclose(logger->logfile)) {
        perror("ERROR: Closing logging file failed");
        exit(EXIT_FAILURE);
      }
    }
    free(logger);
  }
}

#define LEVEL_STR(LEVEL) \
  ((LEVEL) == LOG_ERROR    ? "ERROR" : \
   (LEVEL) == LOG_WARN     ? "WARNING" : \
   (LEVEL) == LOG_SUCCESS  ? "SUCCESS" : \
   (LEVEL) == LOG_PHASE    ? "PHASE" : \
   (LEVEL) == LOG_INFO     ? "INFO" : \
   (LEVEL) == LOG_PARAM    ? "PARAM" : \
   (LEVEL) == LOG_MEM      ? "MEMORY" : \
   (LEVEL) == LOG_TRACE    ? "TRACE" : \
                             "DEBUG")

#define LEVEL_COLOR(LEVEL) \
  ((LEVEL) == LOG_ERROR    ? COLOR_BRIGHT_RED     : \
   (LEVEL) == LOG_WARN     ? COLOR_BRIGHT_MAGENTA : \
   (LEVEL) == LOG_SUCCESS  ? COLOR_BRIGHT_GREEN   : \
   (LEVEL) == LOG_PHASE    ? COLOR_BRIGHT_YELLOW  : \
   (LEVEL) == LOG_INFO     ? COLOR_YELLOW         : \
   (LEVEL) == LOG_PARAM    ? COLOR_BRIGHT_BLUE    : \
   (LEVEL) == LOG_TRACE    ? COLOR_BRIGHT_BLACK   : \
   (LEVEL) == LOG_MEM      ? COLOR_MAGENTA        : \
                             COLOR_CYAN)

/* localtime_r: the loggers can be used by several threads */
#define GET_TIME(S) {\
  time_t now = time(0); \
  struct tm now_tm; \
  if (localtime_r(&now, &now_tm) == NULL || \
      strftime(S, sizeof(S), "%F %T", &now_tm) == 0) \
    S[0] = '\0'; \
}

#define LINE_PREFIX(LOGFILE, LEVEL) {\
  char now_str[100]; \
  GET_TIME(now_str); \
  fprintf(LOGFILE, "# %s | %s | ", now_str, LEVEL_STR(LEVEL)); \
}

#define LINE_DEBUG_PREFIX(LOGFILE, DFILE, DLINE) {\
  fprintf(LOGFILE, "%s:%i | ", DFILE, DLINE); \
}

#define LPRINTF(STREAM, COLORED, PREFIXED, LEVEL, FMT, DBG, DFILE, DLINE) {\
  va_list arglist;\
  if (COLORED)\
    START_COLOR(STREAM, LEVEL_COLOR(LEVEL));\
  if (PREFIXED) \
  {\
    LINE_PREFIX(STREAM, LEVEL);\
    if (DBG) \
      LINE_DEBUG_PREFIX(STREAM, DFILE, DLINE);\
  }\
  va_start(arglist, FMT);\
  vfprintf(STREAM, FMT, arglist);\
  if (COLORED)\
    STOP_COLOR(STREAM);\
  va_end(arglist);\
}

void _fatal_printf(const char *file, int line, const char *fmt, ...) {
  bool colored = isatty(fileno(stderr));
  LPRINTF(stderr, colored, true, LOG_ERROR, fmt, true, file, line);
  exit(EXIT_FAILURE);
}

void _TRACE(const char *file, int line, const char *fmt, ...) {
  bool colored = isatty(fileno(stderr));
  LPRINTF(stderr, colored, true, LOG_TRACE, fmt, true, file, line);
}

void NTRACE(const char *fmt, ...) {
  bool colored = isatty(fileno(stderr));
  LPRINTF(stderr, colored, false, LOG_TRACE, fmt, false, "", 0);
}

/* level of the messages printed (to stderr) without a logger (NULL) */
#define DEFAULT_LOGGER_LEVEL LOG_ERROR

void _lprintf(Logger *logger, LoggingLevel level, const char *file,
     int line, const char *fmt, ...) {
  if (logger) {
    if (level <= logger->level)
      LPRINTF(logger->logfile, logger->colors, true, level, fmt,
              level >= LOG_DEBUG, file, line);
  }
  else if (level <= DEFAULT_LOGGER_LEVEL) {
    bool colored = isatty(fileno(stderr));
    LPRINTF(stderr, colored, true, level, fmt, level >= LOG_DEBUG, file, line);
  }
}

void lnprintf(Logger *logger, LoggingLevel level, const char *fmt, ...) {
  if (logger) {
    if (level <= logger->level)
      LPRINTF(logger->logfile, logger->colors, false, level, fmt,
              false, "", 0);
  }
  else if (level <= DEFAULT_LOGGER_LEVEL) {
      bool colored = isatty(fileno(stderr));
      LPRINTF(stderr, colored, false, level, fmt, false, "", 0);
  }
}

void _showmem(Logger *logger, const char *file, int line) {
  struct rusage r_usage;
  if (getrusage(RUSAGE_SELF, &r_usage) != 0) {
    _lprintf(logger, LOG_MEM, file, line,
              "Memory usage (HWM): not available, %s\n", strerror(errno));
    return;
  }
#ifdef __APPLE__
  /* in bytes on macOS, in kilobytes on Linux */
  r_usage.ru_maxrss /= 1024;
#endif
  _lprintf(logger, LOG_MEM, file, line,
      "Memory usage (HWM): max=%ld kilobytes\n",
      r_usage.ru_maxrss);
}

