/* Best affine projective map from many point pairs. */
/* Last edited on 2024-12-05 10:19:25 by stolfi */ 
   
#ifndef hr2_pmap_aff_from_point_pairs_H
#define hr2_pmap_aff_from_point_pairs_H

#include <r2.h>
#include <stdint.h>
#include <r3.h>
#include <r2x2.h>
#include <r3x3.h>
#include <sign.h>
     
hr2_pmap_t hr2_affmap_from_point_pairs(int32_t np, r2_t p1[], r2_t p2[]);
  /* Returns an affine map {M} (with {M[1,0]=M[2,0]=0}) that best maps
    the points {p1[0..np-1]} to the points {p2[0..np-1]}, in the sense of minimizing the
    mean squared error.
    
    If {np==0}, the reslt is the identity map. If {np==1}, the result is
    a translation. If {np==2}, the result is an Euclidean similarity
    (translation, rotation and scaling). If {np} is 3 or more the result
    is a general affine map. */

#endif
