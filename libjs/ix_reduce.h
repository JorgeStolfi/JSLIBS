#ifndef ix_reduce_H
#define ix_reduce_H

/* Reducing out-of-bounds indices to given range. */
/* Last edited on 2024-11-23 05:52:15 by stolfi */

#include <stdint.h>

#include <ix_types.h>

/* INDEX RANGE REDUCTION */

typedef enum
  { ix_reduce_mode_SINGLE,  /* ... *,*,*,0,1,2,3,4,5,*,*,*,*,*,*,*,*,*,* ... */
    ix_reduce_mode_EXTEND,  /* ... 0,0,0,0,1,2,3,4,5,5,5,5,5,5,5,5,5,5,5 ... */
    ix_reduce_mode_REPEAT,  /* ... 3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3 ... */
    ix_reduce_mode_MIRROR,  /* ... 2,1,0,0,1,2,3,4,5,5,4,3,2,1,0,0,1,2,3 ... */
    ix_reduce_mode_PXMIRR   /* ... 3,2,1,0,1,2,3,4,5,4,3,2,1,0,1,2,3,4,5 ... */
  } ix_reduce_mode_t;
#define ix_reduce_mode_FIRST ix_reduce_mode_SINGLE
#define ix_reduce_mode_LAST  ix_reduce_mode_PXMIRR
  /* A mapping from unrestricted {ix_index_t} values (positive or negative) to 
    integers in some range {0..N-1}.  The comments above illustrate each mapping
    for {N=6} and various indices from {-3} to {+15}.  The code '*' denotes {-1}. */

ix_index_t ix_reduce ( ix_index_t i, ix_size_t N, ix_reduce_mode_t red );
  /* Returns the index {i} reduced to {0..N-1} as specified by {red}.
    The range size {N} must be positive and at most {INT64_MAX}
    (not {UINT64_MAX}). */

#endif
