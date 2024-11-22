/* Binary search on a vector of {int32_t} elements */
/* Last edited on 2024-11-15 19:11:32 by stolfi */

#ifndef binsearch_int32_H
#define binsearch_int32_H

#include <stdint.h>
#include <stdlib.h>

uint64_t binsearch_int32(int32_t y, uint64_t n, int32_t x[]);
  /* Binary search of {y} in {x[0..n-1]}.
    Assumes that the integers {x[0..n-1]} are sorted in non-decreasing order.
    Returns an index {k} such that {x[k-1] < y <= x[k]}, if there is such {k}
    in {1..n-1}. Returns 0 if {y <= x[0]} (or if {n} is zero). 
    Returns {n} if {x[n-1] < y}. */ 

#endif
