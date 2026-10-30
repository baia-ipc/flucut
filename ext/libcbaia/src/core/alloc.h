/* (c) 2017-2021 Giorgio Gonnella */

/* Memory allocation macros*/

#ifndef ALLOC_H
#define ALLOC_H

#include <stdio.h>
#include <stdlib.h>
#include "logger.h"

/* Allocate a buffer of the given size, exit program if failed. */
#define XMALLOC(PTR, BUFSIZE) do {\
  PTR = malloc((BUFSIZE)); \
  if (PTR == NULL) { \
    fatal_printf("Failed allocating %zu bytes\n", (size_t)(BUFSIZE)); \
  }\
  DEBUG_TRACE("Allocated %zu bytes for %s\n",\
         (size_t)(BUFSIZE), #PTR);\
} while (0)

/* Allocate and initialize buffer of the given size, exit program if failed. */
#define XCALLOC(PTR, NITEMS, ITMSIZE) do {\
  PTR = calloc((NITEMS), (ITMSIZE)); \
  if (PTR == NULL) { \
    fatal_printf("Failed allocating %zu bytes\n",\
                 (size_t)(ITMSIZE) * (NITEMS)); \
  }\
  DEBUG_TRACE("Allocated %zu bytes for %s\n",\
        (size_t)(NITEMS) * (ITMSIZE), #PTR);\
} while (0)

/* Reallocate buffer of the given size, exit program if failed. */
#define XREALLOC(PTR, REQSIZE) do { \
  PTR = realloc(PTR, (REQSIZE));\
  if (PTR == NULL) { \
    fatal_printf("Failed allocating %zu bytes\n", (size_t)(REQSIZE)); \
  }\
  DEBUG_TRACE( "Re-allocated %s to %zu bytes\n",\
         #PTR, (size_t)(REQSIZE));\
} while (0)

/* Reallocate a buffer if there is not enough space (to exactly the required
   size, in bytes). Note: PTRPTR is a pointer to the buffer pointer, while
   the other XCHECKSPACE* macros take the buffer pointer itself. */
#define XCHECKSPACE(PTRPTR, SIZEPTR, REQSIZE) do { \
  if ((REQSIZE) > *(SIZEPTR)) {\
    XREALLOC((*(PTRPTR)), (REQSIZE));\
    (*(SIZEPTR)) = (REQSIZE);\
  }\
} while (0)

/* Reallocate a buffer if there is not enough space (to the required size
   plus INC bytes). */
#define XCHECKSPACEINC(PTR, SIZEPTR, REQSIZE, INC) do { \
  if ((REQSIZE) > *(SIZEPTR)) {\
    XREALLOC(PTR, (REQSIZE) + (INC));\
    (*(SIZEPTR)) = (REQSIZE) + (INC);\
  }\
} while (0)

/* Reallocate a buffer if there is not enough space (to the required size
   plus INCPERC percent of it, for amortized constant time appends). */
#define XCHECKSPACEINCPERC(PTR, SIZEPTR, REQSIZE, INCPERC) do { \
  XCHECKSPACEINC(PTR, SIZEPTR, REQSIZE, (((INCPERC)*(REQSIZE))/100));\
} while (0)

/* Verbose version of free() */
#define XFREE(PTR) do { \
  free(PTR);\
  PTR = NULL;\
  DEBUG_TRACE("File %s Line %i | "\
         "Freed memory for %s\n",\
         __FILE__, __LINE__, #PTR);\
} while (0)

#endif
