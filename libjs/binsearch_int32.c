
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>

#include <binsearch_int32.h>

uint64_t binsearch_int32(int32_t y, uint64_t n, int32_t x[])
  { uint64_t i = 0, j = n;
    /* Invariant (A): 
        {0 <= i <= j <= n} /\
        {x[i-1] < y <= x[j]}
        assuming {x[-1] = -oo} and {x[n] = +oo}.
    */
    while (i < j)
      { /* Invariant: (A) and {i < j}. */ 
        uint64_t k = (i + j)/2; /* Cannot be {n}. */
        assert((i <= k) && (k < j));
        if (y <= x[k])
          { j = k; }
        else
          { i = k+1; }
      }
    assert(i == j);
    return i;
  }
