/* (c) 2018 Giorgio Gonnella */

/*
   Implementation:
   - Sorted array S plus an unsorted buffer B
   - Query is by binary search on S plus linear search on B
   - Insertion is in B
   - B is sorted and merged in S when full
*/

#ifndef SET_H
#define SET_H

#include <stdbool.h>
#include <stddef.h>

typedef size_t set_elem_t;

typedef struct Set Set;

/* Create a new set

   <insertion_buffer_size> size of the unsorted buffer used on insertion
                           in number of elements (at least 1, the program
                           is exited with an error message otherwise)

   <bloom_filter_size> if not zero: speedup queries by using a bloom filter
                                    of size 2 ** bloom_filter_size bytes

   <k> if not zero: set is composed of 2**k buckets and the elements are
                    added to/queried from these buckets based on their
                    lowest k bits (k smaller than the number of bits of
                    set_elem_t)
*/

Set *set_new(size_t insertion_buffer_size, size_t bloom_filter_size,
             size_t k);

/* Inserts an element in the set.

   Note: It is not checked if the element is already
   present in the set.
*/
void set_insert_simple(Set *set, set_elem_t key);

/* Inserts an element in the set.

   Note: It is checked if the element is already
   present in the set, in which case it is not inserted.
*/
void set_insert(Set *set, set_elem_t key);

/* Modify the internal state after multiple insertions in order
   to speed up subsequent queries.

   Details: the _insert/_insert_simple inserts an element in an unsorted
   insertion buffer, in which a linear search is applied by queries. Instead,
   the main buffer is sorted and the queries work by binary search. This method
   sorts the insertion buffer and merges it to the main buffer (which is
   otherwise done when the insertion buffer is full).
 */
void set_finalize_insertions(Set *set);

/* Check if an element is present in the set */
bool set_query(const Set *set, set_elem_t key);

/* Delete a set (NULL is accepted) */
void set_delete(Set *set);

/* Unit tests */
void set_test(void);

#endif
