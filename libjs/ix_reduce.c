/* See {ix_reduce.h} */
/* Last edited on 2024-11-23 06:09:57 by stolfi */

#include <stdint.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>

#include <ix_types.h>
#include <ix_reduce.h>

ix_index_t ix_reduce(ix_index_t i, ix_size_t N, ix_reduce_mode_t red)
  {
    demand((N > 0) && (N <= INT64_MAX), "{N} too large");
    if ((i < 0) || (i >= N))
      { /* Reduction is needed. */
        /* Beware that {s % u} with {s} signed and {u} unsigned is {((unsiged)s) % u}. */
        if (red == ix_reduce_mode_SINGLE)
          { i = -1; }
        else if (red == ix_reduce_mode_EXTEND)
          { i = (i < 0 ? 0 : (int64_t)N-1); }
        else if (N == 1)
          { i = 0; }
        else
          { int64_t M = (red == ix_reduce_mode_PXMIRR ? (int64_t)N-1 :(int64_t)N);
            /* Reduce {i} modulo {M}: */
            assert(M >= 1);
            int64_t q = i / M;
            if (i < 0) { q--; }
            i = i - q*M;
            assert((i >= 0) && (i < M));
            
            switch(red)
              { case ix_reduce_mode_REPEAT:  /* ... 3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3 ... */
                  /* Nothing to do: */
                  break;
                case ix_reduce_mode_MIRROR:  /* ... 2,1,0,0,1,2,3,4,5,5,4,3,2,1,0,0,1,2,3 ... */
                case ix_reduce_mode_PXMIRR:  /* ... 3,2,1,0,1,2,3,4,5,4,3,2,1,0,1,2,3,4,5 ... */
                  if ((q & 1) == 1) { i = (int64_t)N - 1 - i ; }
                  break;
                default:
                  demand(FALSE, "invalid reduction code");
              }
          }
      }
    return i;
  }
