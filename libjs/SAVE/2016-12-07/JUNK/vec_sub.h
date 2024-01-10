/* Subsampled vectors. */
/* Last edited on 2008-02-23 09:14:23 by stolfi */

#ifndef vec_sub_H
#define vec_sub_H

#include <stdlib.h>
#include <stdint.h>

#include <bool.h>
#include <ref.h>
#include <vec.h>

/* SUB-VECTOR DESCRIPTORS */

typedef struct vec_sub_t /* Regular sub-sequence of a {vec_t}. */
  { uint32_t nel;   /* Number of elements in sub-vector. */
    void *el;  /* Pointer to element zero of sub-vector. */
    int step;  /* Increment between consecutive elements (in elements). */
  } vec_sub_t;

vec_sub_t vec_as_sub(vec_t *v);
  /* Returns a sub-vector descriptor for the entire vector
    {v} (with {step == 1}). */
    
vec_sub_t vec_sub(vec_sub_t *vs, uint32_t start, int step, uint32_t nel, size_t elsz);
  /* Returns a descriptor {r} for the specified subsequence of elements
    of {vs}, namely {r[i] = vs[start + i*step]}, for {i = 0..r.nel-1}.
    The {step} cannot be zero.
    
    The result size {r.nel} is the largest integer in {0..nel} such
    that all those elements exist. In particular, if {start} is not in
    {0..vs.nel-1}, then {m} is zero. */

#endif
