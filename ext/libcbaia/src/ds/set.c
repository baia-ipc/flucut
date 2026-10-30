/* (c) 2018 Giorgio Gonnella */

#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include "alloc.h"
#include "bloom-filter.h"
#include "set.h"

/*
  Compilation flags:
  SET_USE_SET_BSEARCH   use set_bsearch (default: use stdlib bsearch)
  GG_SET_DEBUG          print the operations and the buffers to stdout
*/

#define NSETS(K) ((size_t)1 << (K))
#define MASK(K) (NSETS(K) - 1)

struct Set
{
  size_t k;
  set_elem_t** mainbuf;
  size_t *mainbufsize;
  size_t *mainbuffill;
  set_elem_t** insbuf;
  size_t *insbufsize;
  size_t *insbuffill;
  bloom_filter_t *bf;
};

#define MAINBUF_INC 16384

/* increment of the size of the main buffer of a bucket: MAINBUF_INC over
   all buckets, but at least 1 (with k >= 15 the shift would give 0) */
static size_t mainbuf_inc(size_t k)
{
  return k < 15 ? (size_t)MAINBUF_INC >> k : 1;
}

Set *set_new(size_t insertion_buffer_size, size_t bloom_filter_size,
             size_t k)
{
  Set *set;
  size_t i;
  assert(k < sizeof(set_elem_t) * CHAR_BIT);
  if (insertion_buffer_size == 0)
    fatal_printf("set_new: the insertion buffer size must be at least 1\n");
  size_t nsets = NSETS(k);
  XMALLOC(set, sizeof *set);
  XMALLOC(set->mainbuf, nsets * sizeof(*set->mainbuf));
  XMALLOC(set->mainbufsize, nsets * sizeof(size_t));
  XMALLOC(set->mainbuffill, nsets * sizeof(size_t));
  XMALLOC(set->insbuf, nsets * sizeof(*set->insbuf));
  XMALLOC(set->insbufsize, nsets * sizeof(size_t));
  XMALLOC(set->insbuffill, nsets * sizeof(size_t));
  for (i = 0; i < nsets; i++) {
    set->mainbufsize[i] = mainbuf_inc(k);
    set->mainbuffill[i] = 0;
    XMALLOC(set->mainbuf[i], set->mainbufsize[i] * sizeof(set_elem_t));
    set->insbufsize[i] = insertion_buffer_size;
    XMALLOC(set->insbuf[i], set->insbufsize[i] * sizeof(set_elem_t));
    set->insbuffill[i] = 0;
  }
  set->k = k;
  set->bf = (bloom_filter_size > 0) ?
    bloom_filter_new(bloom_filter_size) : NULL;
#ifdef GG_SET_DEBUG
  printf("Created set %p\n", (void *)set);
#endif
  return set;
}

static int set_elem_cmp(const void *a_ptr, const void *b_ptr)
{
  set_elem_t a = *(set_elem_t*)a_ptr, b = *(set_elem_t*)b_ptr;
#ifdef GG_SET_DEBUG
  printf("set_elem_cmp; keys %zu and %zu = %d\n", a, b, (a > b) - (a < b));
#endif
  return (a > b) - (a < b);
}

#ifdef GG_SET_DEBUG
static void dbg_show_buffer(const Set *set, size_t i, bool ins, const char *msg)
{
  size_t j, fill;
  set_elem_t *buf;
  const char *bufname;
  if (ins)
  {
    fill = set->insbuffill[i];
    buf = set->insbuf[i];
    bufname = "insbuf";
  }
  else
  {
    fill = set->mainbuffill[i];
    buf = set->mainbuf[i];
    bufname = "mainbuf";
  }
  printf("%s[%zu] %s (fill=%zu): ", bufname, i, msg, fill);
  for (j = 0; j < fill; j++)
  {
    if (j > 0)
      printf(" -- ");
    printf("%zu", buf[j]);
  }
  printf("\n");
}
#endif

static void merge_buffers(Set *set, size_t i)
{
  size_t finalfill, finalpos;
#ifdef GG_SET_DEBUG
  dbg_show_buffer(set, i, false, " before merging");
#endif
  finalfill = set->mainbuffill[i] + set->insbuffill[i];
  assert(set->mainbufsize[i] >= finalfill);
  for (finalpos = finalfill; finalpos > 0; finalpos--)
  {
    assert(finalpos < set->mainbufsize[i]);
    if (set->insbuffill[i] == 0)
    {
      assert(set->mainbuffill[i] == finalpos);
      set->mainbuffill[i] = 0;
      break;
    }
    else if (set->mainbuffill[i] == 0)
    {
      assert(set->insbuffill[i] > 0);
      set->mainbuf[i][finalpos - 1] = set->insbuf[i][set->insbuffill[i] - 1];
      (set->insbuffill[i])--;
    }
    else if (set->insbuf[i][set->insbuffill[i] - 1] >
             set->mainbuf[i][set->mainbuffill[i] - 1])
    {
      set->mainbuf[i][finalpos - 1] = set->insbuf[i][set->insbuffill[i] - 1];
      (set->insbuffill[i])--;
    }
    else
    {
      set->mainbuf[i][finalpos - 1] = set->mainbuf[i][set->mainbuffill[i] - 1];
      (set->mainbuffill[i])--;
    }
  }
  assert(set->insbuffill[i] == 0);
  assert(set->mainbuffill[i] == 0);
  set->mainbuffill[i] = finalfill;
#ifdef GG_SET_DEBUG
  dbg_show_buffer(set, i, false, " after merging");
#endif
}

static void set_finalize_insertions_i(Set *set, size_t i)
{
#ifdef GG_SET_DEBUG
  printf("sorting insbuf %zu\n", i);
  dbg_show_buffer(set, i, true, " before sorting");
#endif
  qsort(set->insbuf[i], set->insbuffill[i],
      sizeof(set_elem_t), set_elem_cmp);
#ifdef GG_SET_DEBUG
  dbg_show_buffer(set, i, true, " after sorting");
#endif
  if (set->mainbuffill[i] + set->insbuffill[i] >= set->mainbufsize[i])
  {
    size_t newsize = set->mainbufsize[i] + set->insbuffill[i] +
                     mainbuf_inc(set->k);
#ifdef GG_SET_DEBUG
    printf("realloc set %p mainbuf; prev_size=%zu, new_size = %zu\n",
        (void *)set, set->mainbufsize[i], newsize);
#endif
    set->mainbufsize[i] = newsize;
    XREALLOC(set->mainbuf[i], newsize * sizeof *set->mainbuf[i]);
  }
  merge_buffers(set, i);
}

void set_finalize_insertions(Set *set) {
  size_t i, nsets = NSETS(set->k);
  for (i = 0; i < nsets; i++)
    set_finalize_insertions_i(set, i);
}

void set_insert_simple(Set *set, set_elem_t key) {
  if (set->bf != NULL)
    bloom_filter_insert(set->bf, key);
  size_t i = key & MASK(set->k);
  key = key >> set->k;
#ifdef GG_SET_DEBUG
  printf("\nInsert key %zu in set %p i=%zu\n", key, (void *)set, i);
  printf("Is=%zu;If=%zu;Ms=%zu;Mf=%zu\n",
      set->insbufsize[i], set->insbuffill[i],
      set->mainbufsize[i], set->mainbuffill[i]);
#endif
  assert(set->insbuffill[i] < set->insbufsize[i]);
  set->insbuf[i][set->insbuffill[i]] = key;
  set->insbuffill[i]++;
  if (set->insbuffill[i] == set->insbufsize[i])
    set_finalize_insertions_i(set, i);
}

void set_insert(Set *set, set_elem_t key)
{
  if (!set_query(set, key))
    set_insert_simple(set, key);
}

#ifdef SET_USE_SET_BSEARCH
static bool set_bsearch(const Set *set, size_t i, set_elem_t key)
{
  if (set->mainbuffill[i] > 0) {
    size_t lower = 0, upper = set->mainbuffill[i]-1, mid;
    set_elem_t midvalue;
    while (lower <= upper) {
      mid = lower + ((upper-lower)>>1);
      assert(mid < set->mainbuffill[i]);
      midvalue = set->mainbuf[i][mid];
      if (key < midvalue) {
        if (mid > 0)
          upper = mid-1;
        else
          return false;
      }
      else if (key > midvalue)
        lower = mid+1;
      else {
        return true;
      }
    }
  }
  return false;
}
#endif

bool set_query(const Set *set, set_elem_t key)
{
  if (set->bf != NULL && !bloom_filter_query(set->bf, key))
    return false;
  size_t i = key & MASK(set->k);
  key = key >> set->k;
#ifdef GG_SET_DEBUG
  printf("\nSearch key: %zu in set %p (i=%zu)\n", key, (void *)set, i);
  printf("Is=%zu;If=%zu;Ms=%zu;Mf=%zu\n",
      set->insbufsize[i], set->insbuffill[i],
      set->mainbufsize[i], set->mainbuffill[i]);
  dbg_show_buffer(set, i, false, "");
  dbg_show_buffer(set, i, true, "");
#endif
#ifdef SET_USE_SET_BSEARCH
  if (set_bsearch(set, i, key))
#else
  if (bsearch(&key, set->mainbuf[i], set->mainbuffill[i],
            sizeof key, set_elem_cmp) != NULL)
#endif
  {
#ifdef GG_SET_DEBUG
    printf("key found in main buffer %zu\n", i);
#endif
    return true;
  }
  else
  {
    size_t j;
    for (j = 0; j < set->insbuffill[i]; j++)
    {
      if (set->insbuf[i][j] == key)
      {
#ifdef GG_SET_DEBUG
        printf("key found in ins buffer %zu\n", i);
#endif
        return true;
      }
    }
#ifdef GG_SET_DEBUG
    printf("key not found\n");
#endif
    return false;
  }
}

void set_delete(Set *set)
{
  if (set == NULL)
    return;
  size_t i, nsets = NSETS(set->k);
  for (i = 0; i < nsets; i++) {
    free(set->mainbuf[i]);
    free(set->insbuf[i]);
  }
  free(set->mainbuf);
  free(set->mainbufsize);
  free(set->mainbuffill);
  free(set->insbuf);
  free(set->insbufsize);
  free(set->insbuffill);
  bloom_filter_delete(set->bf);
  free(set);
}

#include "test.h"

static size_t set_mainbuffill(Set *set) {
  size_t i, nsets = NSETS(set->k), result = 0;
  for (i = 0; i < nsets; i++)
    result += set->mainbuffill[i];
  return result;
}

static size_t set_insbuffill(Set *set) {
  size_t i, nsets = NSETS(set->k), result = 0;
  for (i = 0; i < nsets; i++)
    result += set->insbuffill[i];
  return result;
}

void set_test(void)
{
  Set *set = set_new(3, 0, 0);
  test_assert(set_mainbuffill(set) == 0);
  test_assert(set_insbuffill(set) == 0);
  test_assert(!set_query(set, 1));
  test_assert(!set_query(set, 2));
  test_assert(!set_query(set, 3));
  test_assert(!set_query(set, 4));
  test_assert(!set_query(set, 5));
  test_assert(!set_query(set, 6));
  test_assert(!set_query(set, 7));
  set_insert(set, 3);
  test_assert(set_mainbuffill(set) == 0);
  test_assert(set_insbuffill(set) == 1);
  test_assert(!set_query(set, 1));
  test_assert(!set_query(set, 2));
  test_assert(set_query(set, 3));
  test_assert(!set_query(set, 4));
  test_assert(!set_query(set, 5));
  test_assert(!set_query(set, 6));
  test_assert(!set_query(set, 7));
  set_insert(set, 1);
  test_assert(set_mainbuffill(set) == 0);
  test_assert(set_insbuffill(set) == 2);
  test_assert(set_query(set, 1));
  test_assert(!set_query(set, 2));
  test_assert(set_query(set, 3));
  test_assert(!set_query(set, 4));
  test_assert(!set_query(set, 5));
  test_assert(!set_query(set, 6));
  test_assert(!set_query(set, 7));
  set_insert(set, 4);
  test_assert(set_mainbuffill(set) == 3);
  test_assert(set_insbuffill(set) == 0);
  test_assert(set_query(set, 1));
  test_assert(!set_query(set, 2));
  test_assert(set_query(set, 3));
  test_assert(set_query(set, 4));
  test_assert(!set_query(set, 5));
  test_assert(!set_query(set, 6));
  test_assert(!set_query(set, 7));
  set_insert(set, 2);
  test_assert(set_mainbuffill(set) == 3);
  test_assert(set_insbuffill(set) == 1);
  test_assert(set_query(set, 1));
  test_assert(set_query(set, 2));
  test_assert(set_query(set, 3));
  test_assert(set_query(set, 4));
  test_assert(!set_query(set, 5));
  test_assert(!set_query(set, 6));
  test_assert(!set_query(set, 7));
  set_insert(set, 6);
  test_assert(set_mainbuffill(set) == 3);
  test_assert(set_insbuffill(set) == 2);
  test_assert(set_query(set, 1));
  test_assert(set_query(set, 2));
  test_assert(set_query(set, 3));
  test_assert(set_query(set, 4));
  test_assert(!set_query(set, 5));
  test_assert(set_query(set, 6));
  test_assert(!set_query(set, 7));
  set_insert(set, 5);
  test_assert(set_mainbuffill(set) == 6);
  test_assert(set_insbuffill(set) == 0);
  test_assert(set_query(set, 1));
  test_assert(set_query(set, 2));
  test_assert(set_query(set, 3));
  test_assert(set_query(set, 4));
  test_assert(set_query(set, 5));
  test_assert(set_query(set, 6));
  test_assert(!set_query(set, 7));
  set_insert(set, 7);
  test_assert(set_mainbuffill(set) == 6);
  test_assert(set_insbuffill(set) == 1);
  test_assert(set_query(set, 1));
  test_assert(set_query(set, 2));
  test_assert(set_query(set, 3));
  test_assert(set_query(set, 4));
  test_assert(set_query(set, 5));
  test_assert(set_query(set, 6));
  test_assert(set_query(set, 7));
  set_delete(set);
  /* with bloom filter and buckets */
  set = set_new(4, 10, 2);
  for (set_elem_t key = 0; key < 1000; key += 3)
    set_insert(set, key);
  set_finalize_insertions(set);
  for (set_elem_t key = 0; key < 1000; key++)
    test_assert(set_query(set, key) == (key % 3 == 0));
  set_delete(set);
  /* many buckets: the initial main buffers are not empty */
  set = set_new(2, 0, 16);
  for (set_elem_t key = 0; key < 200000; key += 7)
    set_insert(set, key);
  set_finalize_insertions(set);
  for (set_elem_t key = 0; key < 200000; key++)
    test_assert(set_query(set, key) == (key % 7 == 0));
  set_delete(set);
}
