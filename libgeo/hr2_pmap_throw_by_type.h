#ifndef hr2_pmap_throw_by_type_H
#define hr2_pmap_throw_by_type_H

/* Generate random projective maps of a given type and handedness. */
/* Last edited on 2024-12-05 10:27:10 by stolfi */ 

#include <stdint.h>

#include <bool.h>
#include <hr2_pmap.h>
#include <hr2.h>

hr2_pmap_t hr2_pmap_throw_by_type(hr2_pmap_type_t type, sign_t sgn);
  /* Returns a random projective map with given {type}.  
    The map will have the handedness {sgn}, which must be {-1} or {+1}.
    The two matrices {M.dir} and {M.inv} will be normalized
    to unit Euclidean norm. */

#endif
