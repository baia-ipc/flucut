/* (c) 2026 Giorgio Gonnella */

#include <assert.h>
#include <stdlib.h>
#include "alloc.h"
#include "heap.h"
#include "test.h"

#define HEAP_MIN_CAPACITY 16

struct heap_t {
  heap_cmp_t cmp;
  void **elements;   /* elements[0] is the root */
  size_t len, cap;
};

heap_t* heap_new(heap_cmp_t cmp, size_t capacity) {
  assert(cmp != NULL);
  heap_t *self;
  XMALLOC(self, sizeof(heap_t));
  self->cmp = cmp;
  self->len = 0;
  self->cap = capacity < HEAP_MIN_CAPACITY ? HEAP_MIN_CAPACITY : capacity;
  XMALLOC(self->elements, sizeof(void *) * self->cap);
  return self;
}

void heap_delete(heap_t *self) {
  if (self == NULL)
    return;
  free(self->elements);
  free(self);
}

size_t heap_size(const heap_t *self) {
  assert(self != NULL);
  return self->len;
}

bool heap_is_empty(const heap_t *self) {
  assert(self != NULL);
  return self->len == 0;
}

void heap_push(heap_t *self, void *element) {
  assert(self != NULL);
  if (self->len == self->cap) {
    self->cap *= 2;
    XREALLOC(self->elements, sizeof(void *) * self->cap);
  }
  /* sift up: move the parents down until the place is found */
  size_t i = self->len++;
  while (i > 0) {
    size_t parent = (i - 1) / 2;
    if (self->cmp(self->elements[parent], element) <= 0)
      break;
    self->elements[i] = self->elements[parent];
    i = parent;
  }
  self->elements[i] = element;
}

void *heap_peek(const heap_t *self) {
  assert(self != NULL);
  assert(self->len > 0);
  return self->elements[0];
}

void *heap_pop(heap_t *self) {
  assert(self != NULL);
  assert(self->len > 0);
  void *min = self->elements[0];
  void *last = self->elements[--self->len];
  /* sift down: move the smaller children up until the place of the last
     element is found */
  size_t i = 0;
  while (true) {
    size_t child = 2 * i + 1;
    if (child >= self->len)
      break;
    if (child + 1 < self->len &&
        self->cmp(self->elements[child + 1], self->elements[child]) < 0)
      child++;
    if (self->cmp(last, self->elements[child]) <= 0)
      break;
    self->elements[i] = self->elements[child];
    i = child;
  }
  if (self->len > 0)
    self->elements[i] = last;
  return min;
}

static int cmp_ints(const void *a, const void *b) {
  int x = *(const int *)a, y = *(const int *)b;
  return (x > y) - (x < y);
}

void heap_test(void) {
  int values[] = {5, 3, 9, 1, 7, 3};
  heap_t *heap = heap_new(cmp_ints, 2);
  test_assert(heap_is_empty(heap));
  for (size_t i = 0; i < 6; i++)
    heap_push(heap, values + i);
  test_eq_size(heap_size(heap), 6);
  test_eq_int(*(int *)heap_peek(heap), 1);
  int expected[] = {1, 3, 3, 5, 7, 9};
  for (size_t i = 0; i < 6; i++)
    test_eq_int(*(int *)heap_pop(heap), expected[i]);
  test_assert(heap_is_empty(heap));
  /* many elements, in pseudo-random order (growing the heap) */
  int many[1000];
  for (int i = 0; i < 1000; i++) {
    many[i] = (i * 7919) % 1009;
    heap_push(heap, many + i);
  }
  int prev = -1;
  while (!heap_is_empty(heap)) {
    int v = *(int *)heap_pop(heap);
    test_assert(v >= prev);
    prev = v;
  }
  heap_delete(heap);
}
