#ifndef hr2_pmap_opt_generic_H
#define hr2_pmap_opt_generic_H

/* Last edited on 2023-10-20 18:44:57 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <hr2.h>

void hr2_pmap_opt_generic_encode(hr2_pmap_t *M, r3x3_t *R, double y[]);
  /* Converts a generic projective map {M} to up to eight parameters {y[0..nv-1]}.
  
    More precisely, let {nv} be the number of non-zero elements in the
    3x3 matrix {R}, and {I} be the 3x3 identity matrix. The vector {y}
    must have at least {nv} elements. Sets {y[0..nv-1]} to {(M[i][j] -
    I[i][j])/R[i][j]} for every {i,j} in {0..2} such that {R[i][j]} is
    nonzero, scanning by rows.
    
    Fails if {R} has any negative element, or more than 8 non-zero ones. */

void hr2_pmap_opt_generic_decode(double y[], r3x3_t *R, hr2_pmap_t *M);
  /* Converts up to 8 parameters {y[0..nv-1]} to a generic map {M}, 
    as the inverse of {hr2_pmap_opt_generic_encode}.
  
    More precisely, let {nv} be the number of non-zero elements in the
    3x3 matrix {R}, and {I} be the 3x3 identity matrix. The vector {y}
    must have at least {nv} elements. Sets {M.dir to the identity, then
    adds {y[k]*R[i][j]} to {M.dir[i][j]} for every {i,j} in {0..2} such
    that {R[i][j]} is nonzero, scanning by rows. Finally sets {M.inv} to
    be the inverse of {M.dir}. */

#endif
