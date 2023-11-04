#ifndef hr2_pmap_opt_affine_H
#define hr2_pmap_opt_affine_H

/* Last edited on 2023-10-20 18:26:08 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <hr2.h>

void hr2_pmap_opt_affine_encode(hr2_pmap_t *M, double rad[], double y[]);
  /* Converts a affine map {M} to six parameters {y[0..5]}.
  
    More precisely, assumes that {M} is a affine map with {M.dir[0][0] =
    1}, and that {y,rad} have at least 6 elements each. Let {D} be the
    3x3 matrix that is the difference of {M.dir} from the identity. Sets
    {y[k]} to {D[i][j]/rad[k]} for {i} in {0..2} and {j} in {1..2}, row
    by row. */

void hr2_pmap_opt_affine_decode(double y[], double rad[], hr2_pmap_t *M);
  /* Converts two parameters {y[0..5]} to a affine map {M}, 
    as the inverse of {hr2_pmap_opt_affine_encode}.
  
    More precisely, assumes that {y,rad} have at least 6 elements each,
    and {M} is a affine map with {M.dir[0][0]=M.inv[0][0]=1}. Let {D} be
    the 3x3 matrix with first column zero and the the other elements,
    row by row, equal to {y[k]*rad[k]}, for {k} in {0..5}. Sets {M.dir}
    to the identity plus {D}, then sets {M.inv} to be the inverse of
    {M.dir}. */

#endif
