/* (c) 2026 Giorgio Gonnella */

/*
   Progress display in the style of the Python tqdm package, e.g.:

     Processing reads: 123456it [00:05, 24691.20it/s]
     Processing reads:  45%|█████▍      | 180000/400000 [00:07<00:09, ...]

   The first form is used if the total is unknown (0). The display is
   refreshed at most every PROGRESS_MIN_INTERVAL seconds, by overwriting
   the current line (carriage return).
*/

#ifndef PROGRESS_H
#define PROGRESS_H

#include <stddef.h>
#include <stdio.h>

#define PROGRESS_MIN_INTERVAL 0.1

typedef struct progress_t progress_t;

/*
   Create a progress display on <file> (e.g. stderr), with the given
   <description> (copied) and <total> (0 if unknown); it is shown
   immediately
*/
progress_t* progress_new(FILE *file, const char *description, size_t total);

/* Add <n> to the progress (the display is refreshed, if needed) */
void progress_update(progress_t *self, size_t n);

/* Show the final state and end the line */
void progress_finish(progress_t *self);

/* Current count */
size_t progress_count(const progress_t *self);

/* Delete the progress display (NULL is accepted) */
void progress_delete(progress_t *self);

/* Unit tests */
void progress_test(void);

#endif
