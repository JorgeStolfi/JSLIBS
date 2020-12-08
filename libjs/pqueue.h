/* pqueue - priority queue of integers with real values. */
/* Last edited on 2017-01-02 21:51:59 by jstolfi */

#ifndef pqueue_H
#define pqueue_H

#define _GNU_SOURCE
#include <stdint.h>
#include <bool.h>

typedef uint32_t pqueue_item_t; 
  /* Items in queue are unsigned integers. */

typedef double pqueue_value_t;
  /* The value of an item in the queue. */

typedef struct pqueue_rep_t pqueue_t;
  /* At any quiescent moment, a priority queue {Q} contains {Q.n}
    distinct /items/, which are non-negative integers. Each item has
    an associated double-precision /value/.
   
    The implementation is designed to allow arbitrary sequences of
    lookups, insertions, deletions, value
    updates, and finding the first item in value order, at low
    amortized cost per operation.
    
    The queue's internal tables are automatically (re)allocated as
    needed when items are inserted. These expansions have only {O(1)}
    cost, amortized per operation. Still, the client may want to avoid
    this overhead by pre-allocating tables of adequate size with
    {pqueue_realloc}.
    
    A {pqueue_t} uses {O(m)} space, where {m} is the largest item that
    was stored into the queue since the last time it was empty. The
    client must not modify any field of {Q} except through the
    procedures in this module. All times below are worst-case; {n}
    means the number of elements currently in the queue. */

#define pqueue_NMAX (1073741824)
  /* Max num of simultaneous items in queue (2^30). */

#define pqueue_ITEM_MAX (1073741823)
  /* Max item (2^30-1). */

/* CREATION */

pqueue_t *pqueue_new(void);
  /* Returns a pointer to a newly created queue structure, 
    intially empty and with increasing (+1) order.  */

/* QUERIES (WITHOUT SIDE EFFECTS) */

typedef uint32_t pqueue_count_t;
  /* A count of item in the queue. */

pqueue_count_t pqueue_count(pqueue_t *Q);
  /* The number of items currently in {Q}. Time: {O(1)}.  */

int pqueue_order(pqueue_t *Q);
  /* The ordering of items in {Q}: +1 means INCREASING (default), 
    -1 means DECREASING. Time: {O(1)}.  */

pqueue_item_t pqueue_head(pqueue_t *Q);
  /* The first item in {Q}, in value order. Namely the item with
    minimum value for order +1, or the maximum value for order = -1.
    Does not change {Q}. Time: {O(1)}.  */

bool_t pqueue_has(pqueue_t *Q, pqueue_item_t z);
  /* TRUE iff item {z} is in {Q}. Time: {O(1)}. */

pqueue_value_t pqueue_value(pqueue_t *Q, pqueue_item_t z);
  /* Current value of item {z} in {Q}; undefined if {z} is not in {Q}.
    Time: {O(1)}. */

/* MODIFYING */

void pqueue_insert(pqueue_t *Q, pqueue_item_t z, pqueue_value_t v);
  /* Inserts the item {z} in {Q}, with value {v}, rearranging the
    contents as needed. Undefined if {z} is already in {Q}.
    
    Time: {O(log(n))} worst-case, if no storage re-allocation is
    needed; otherwise {O(n)} worst-case, but {O(log(n))} amortized. */

void pqueue_delete(pqueue_t *Q, pqueue_item_t z);
  /* Removes the item {z} from {Q}. Undefined if {z}
    is not in {Q}. Time: {O(log(n))}. */

void pqueue_set_value(pqueue_t *Q, pqueue_item_t z, pqueue_value_t v);
  /* Sets the value of item {z} to {v}. Time: {O(log(n))}. */

void pqueue_set_order(pqueue_t *Q, int order);
  /* Specifies the ordering of elements for {pqueue_head}.
    Time: {O(n log n)}. */

void pqueue_reset(pqueue_t *Q);
  /* Deletes all elements.  Time: {O(1)}. Does not free any internal 
    storage. */

void pqueue_free(pqueue_t *Q);
  /* Equivalent to {pqueue_reset(Q); pqueue_realloc(Q,0,0); free(Q)}. */

typedef pqueue_value_t pqueue_value_function_t(pqueue_item_t z, pqueue_value_t v);
  /* A client procedure to recompute an item's value. */

void pqueue_set_all_values(pqueue_t *Q, pqueue_value_function_t *f);
  /* Calls {pqueue_set_value(Q, z, v')} for every item {z} in {Q},
    where {v' = f(z, pqueue_value(Q,z))}.  Time: {O(n*log(n))}. */

/* ALLOCATION HINTS */

void pqueue_realloc(pqueue_t *Q, pqueue_count_t nmax, pqueue_item_t zlim);
  /* The queue {Q} must be empty; (re-)allocates its internal storage
    areas so that at least {nmax} items in the range {0..zlim-1} can be
    inserted before reallocation is needed again. 
    
    If the queue currently has space allocated for twice these amounts
    or more, the storage is freed and reallocated with the proper size.
    In particular, {pqueue_realloc(Q,0,0)} frees all internal
    storage (but not {*Q} itself). Time: {O(1)}. */

/* HANDY TOOLS */

int pqueue_dblcmp(pqueue_value_t x, pqueue_value_t y);
  /* Returns -1, 0, or +1 depending on whether the value {x} is less than,
    equal to, or greater than {y}. */

/* ITEM POSITIONS

    At any moment, each item in the queue also has a distinct
    /position/, an integer {p} in the range {0..Q.n-1}.
    
    Note that the position is NOT the item's rank in the order. The
    position of an item {z} in {Q} may be affected by insertions,
    deletions, or value updates, of {z} or of any other item.
 
    The queue is partially sorted by the value fields, either in
    increasing or decreasing order, according to the binary heap rule;
    namely, the first item in value order is at position 0, and the
    item at any position {p > 0} follows item at position {(p-1)/2} in
    that order. */

typedef uint32_t pqueue_position_t;

pqueue_position_t pqueue_position(pqueue_t *Q, pqueue_item_t z);
  /* Returns the index {ix} in {0..Q.n-1} such that {Q->el[ix] == z};
     or {Q.n} if no such {ix} exists. Time: {O(1)}. */

pqueue_item_t pqueue_item(pqueue_t *Q, pqueue_position_t p);
  /* The item currently at position {p} of {Q}; 
    undefined if {p} is not in {0..Q.n-1}. Time: {O(1)}. */

#endif
