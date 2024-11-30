#ifndef ix_types_H
#define ix_types_H

/* data type defintions for indexing operations in {ix.h}. */
/* Last edited on 2024-11-23 05:48:12 by stolfi */

#include <stdint.h>

typedef uint8_t ix_dim_t;
  /* Type for the number of indices {d}. */

typedef uint8_t ix_axis_t;
  /* Type for axis number, or index of an index, usually from 0 to {d-1}. */

typedef uint64_t ix_size_t;
  /* Type for the size {sz[i]} of an array along some axis. */

typedef int64_t ix_index_t;
  /* Type for individual array element indices {ix[i]}.
    They may be negative on occasion (e.g. in reverse 'for' loops).  */

typedef int64_t ix_step_t; 
  /* Type for position steps {st[i]} along any axis. It is signed
    to allow for array flipping. */

typedef uint64_t ix_pos_t; 
  /* Type for element positions.  It must be always in the range {0..N-1}
    where {N} is the number of elements spanned by the original array. */

#define ix_pos_NONE (UINT64_MAX)
  /* An {ix_pos_t} value that means "no such element". */

typedef uint64_t ix_count_t;
  /* A count of the number of elements in an array. */

#endif
