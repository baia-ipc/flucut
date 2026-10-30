/* (c) 2021 Giorgio Gonnella */

/* Logging library */

#ifndef LOGGER_H
#define LOGGER_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct Logger Logger;

typedef enum LoggingLevel {
  LOG_ERROR   = 0, /* always displayed, error = leading to program exit */
  LOG_WARN    = 1, /* displayed if not -q, warning = not leading to program exit */
  LOG_SUCCESS = 2, /* displayed if not -q, success of an operation */
  LOG_PHASE   = 3, /* displayed if -v, start of an operation */
  LOG_INFO    = 4, /* displayed if -v, general verbose information */
  LOG_PARAM   = 5, /* displayed if -v, input parameters */
  LOG_MEM     = 6, /* displayed if -D, memory high water mark */
  LOG_DEBUG   = 7, /* displayed if -D, debug information */
  LOG_TRACE   = 8  /* used by the [DEBUG_]TRACE/[DEBUG_]NTRACE functions */
} LoggingLevel;

Logger* new_logger(FILE *logfile, LoggingLevel level);

void delete_logger(Logger* logger);

/* printf-like function for the logger;
   - message only displayed if the given level is <= the level of the logger;
   - outputs an header and the message */
#define lprintf(LOGGER, LEVEL, FMT, ...) \
  _lprintf(LOGGER, LEVEL, __FILE__, __LINE__, FMT, ##__VA_ARGS__)

/* printf-like function for the logger;
   - message only displayed if the given level is <= the level of the logger;
   - outputs a message, but no header
   - used after lprintf, e.g. when iterating over some values to display,
     which require multiple prints for a single line  */
void lnprintf(Logger *logger, LoggingLevel level, const char *fmt, ...)
     __attribute__ ((format (printf,3,4)));

/* printf-like macro, which checks a condition: if not holding,
   prints an error message and exits the program with code EXIT_FAILURE */
#define LCHECK(CONDITION, LOGGER, FMT, ...) \
   _LCHECK(CONDITION, LOGGER, __FILE__, __LINE__, FMT, ##__VA_ARGS__)

/* printf-like function, printing a last error message to stderr (formatted
   as by a logger) and exiting the program with code EXIT_FAILURE;
   - used when a logger is not available e.g. used by macros in alloc.h */
#define fatal_printf(FMT, ...) \
  _fatal_printf(__FILE__, __LINE__, FMT, ##__VA_ARGS__)

/* show the amount of memory used by the program until now, MEM level */
#define showmem(logger) \
  _showmem(logger, __FILE__, __LINE__)

/* ---- debugging functionality ---- */

/* TRACE, NTRACE, DEBUG_TRACE and DEBUG_NTRACE are upper case, also where
   they are functions (NTRACE), because they behave like the macros (TRACE,
   DEBUG_TRACE) whose lines they complete */

/* - to print something during debugging, e.g. in a
     #ifdef GG_DEBUG_MODULE;
     uses settings for level LOG_TRACE */
#define TRACE(FMT, ...) \
  _TRACE(__FILE__, __LINE__, FMT, ##__VA_ARGS__)

/* like TRACE, but without line header, to complete the line */
void NTRACE(const char *fmt, ...)
     __attribute__ ((format (printf,1,2)));

#ifdef SHOW_TRACES
/* to print something during debugging, when compiling with flag SHOW_TRACES;
   as TRACE */
#define DEBUG_TRACE(FMT, ...) \
  _TRACE(__FILE__, __LINE__, FMT, ##__VA_ARGS__)

/* like DEBUG_TRACE, but without line header, to complete the line; as
   NTRACE */
#define DEBUG_NTRACE(FMT, ...) NTRACE(FMT, ##__VA_ARGS__)
#else
/* without SHOW_TRACES, DEBUG_[N]TRACE are deleted from the code */
#define DEBUG_TRACE(fmt, ...) ((void)0)
#define DEBUG_NTRACE(fmt, ...) ((void)0)
#endif

/* ---- versions requiring __FILE__ and __LINE__ ---- */

/* the _ versions of the functions require to pass the current file and
   line, which is then automated by macros with the same name wo the _ */

void _lprintf(Logger *logger, LoggingLevel level, const char *file,
     int line, const char *fmt, ...)
     __attribute__ ((format (printf,5,6)));

#define _LCHECK(CONDITION, LOGGER, F, L, FMT, ...) do { \
  if (!(CONDITION)) {\
    _lprintf(LOGGER, LOG_ERROR, F, L, FMT, ##__VA_ARGS__);\
    delete_logger(LOGGER);\
    exit(EXIT_FAILURE);\
  }\
} while (0)

void _fatal_printf(const char *file, int line, const char *fmt, ...)
     __attribute__ ((format (printf,3,4), noreturn));

void _showmem(Logger *logger, const char *file, int line);

void _TRACE(const char *file, int line, const char *fmt, ...)
     __attribute__ ((format (printf,3,4)));

#endif
