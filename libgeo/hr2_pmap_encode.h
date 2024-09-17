#ifndef hr2_pmap_encode_H
#define hr2_pmap_encode_H

/* Tools for minimal encoding of projective maps. */
/* Last edited on 2024-09-17 16:29:27 by stolfi */ 

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <hr2_pmap.h>
#include <hr2.h>
 
int32_t hr2_pmap_encode_num_parameters(hr2_pmap_type_t type);
  /* Number of parameters used to encode a projective map of the specified
    {type}. Namely, 0 for identity, 2 for translation, 3 for congruence,
    4 for similarity, 6 for general affine map, and 9 for general
    projective transformation.  

    Note that in the last case the result is 1 more that the degrees of
    freedom in the map. Therefore, when optimizing a function of the
    map, an extra constraint or penalty term must be added to keep the
    optimum unique. */

void hr2_pmap_encode(hr2_pmap_t *M, hr2_pmap_type_t type, int32_t ny, double y[]);
  /* Converts the map {*M} into the parameters {y[0..ny-1]}. Assumes
    that {M} is of the given {type}; the result is undefined otherwise.
    The number {ny} must be {hr2_pmap_encode_num_parameters(type)}. The
    handedness of {M} is lost. */

void hr2_pmap_decode(int32_t ny, double y[], hr2_pmap_type_t type, sign_t sgn, hr2_pmap_t *M);
  /* Converts the encoded parameters {y[0..ny-1]} into the projective
    map {*M}, with handedness {sgn}. The number {ny} must be
    {hr2_pmap_encode_num_parameters(type)}.
    
    The {sgn} parameter must be {-1} or {+1}, and specifies the
    handedness of the resultimg map. If {sgn} is {-1}, the result is the
    composition {S N} of the XY-swap map with the map {N} obtained with
    {sgn=+1}. */

#endif
