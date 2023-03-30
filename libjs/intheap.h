/* intheap - sorted heap of integers. */
/* Last edited on 2023-03-18 11:27:12 by stolfi */

#ifndef intheap_H
#define intheap_H

#define _GNU_SOURCE
#include <stdint.h>

/* 
  HEAPS OF INTEGERS
  
  The following procedures maintain a heap of integers, stored as an
  integer vector {h[0..n-1]}. The root of the heap is {h[0]}, and the
  children of an item {h[i]} are {h[2*i+1]} and {h[2*i+2]}, iff those
  elements exist. 
  
  The elements are ordered so that {sgn*cmp(h[i],h[j]) <= 0} for every
  item {h[i]} and each child {h[j]} of {h[i]}. The values returned by
  {cmp} should be consistent with an order relation (anti-symmetric,
  reflexive, and transitive). The {sgn} parameter is a client
  convenience feature, and is usually {+1} or {-1}. Neither {cmp} nor
  {sgn} should be changed while the heap contains more than one
  item. */

void ihp_heap_insert(int32_t *h, int32_t *n, int32_t v, int32_t cmp(int32_t x, int32_t y), int32_t sgn);
  /* Inserts the integer {v} in the heap {h[0..n-1]}, rearranging 
   its contents as needed.  Also increments {n}. 
   Takes {O(log(n))} time in the worst case. */

int32_t ihp_heap_pop(int32_t *h, int32_t *n, int32_t cmp(int32_t x, int32_t y), int32_t sgn);
  /* Removes the root item {h[0]} from {h}, rearranging its contents
    as needed. Returns the removed item. Also decrements {n}.  
    Takes {O(log(n))} time in the worst case.  */

#endif
