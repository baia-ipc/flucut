/* (c) 2026 Giorgio Gonnella */

/*
   Binary min-heap of pointers, ordered by a comparison function (the
   smallest element is extracted first), e.g. for priority queues or for
   merging sorted sequences.

   Implementation: array, growing when needed.
*/

#ifndef HEAP_H
#define HEAP_H

#include <stdbool.h>
#include <stddef.h>

/* comparison function, as for qsort(): negative if a < b, 0 if a = b,
   positive if a > b; it is called with the pointers stored in the heap */
typedef int (*heap_cmp_t)(const void *a, const void *b);

typedef struct heap_t heap_t;

/* Create an empty heap (<capacity>: initial capacity, it grows if needed) */
heap_t* heap_new(heap_cmp_t cmp, size_t capacity);

/* Add an element */
void heap_push(heap_t *self, void *element);

/* Remove and return the smallest element (the heap must not be empty) */
void *heap_pop(heap_t *self);

/* Return the smallest element (the heap must not be empty) */
void *heap_peek(const heap_t *self);

/* Number of elements */
size_t heap_size(const heap_t *self);

/* Is the heap empty? */
bool heap_is_empty(const heap_t *self);

/* Delete the heap (not the elements; NULL is accepted) */
void heap_delete(heap_t *self);

/* Unit tests */
void heap_test(void);

#endif
