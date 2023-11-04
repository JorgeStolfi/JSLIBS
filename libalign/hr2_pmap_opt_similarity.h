#ifndef hr2_pmap_opt_similarity_H
#define hr2_pmap_opt_similarity_H

/* Last edited on 2023-10-20 16:30:57 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <hr2.h>

void hr2_pmap_opt_similarity_encode(hr2_pmap_t *M, double rad[], double y[]);
  /* Converts a similarity map {M} to four parameters {y[0..3]}.
  
    More precisely, assumes that {M} is a similarity map
    (combination of translation, rotation, and uniform scaling)
    with {M.dir[0][0]=1}, and that {y,rad} have at least 4 elements each.
    Sets {y[0]} and {y[1]} to {M.dir[0][1]/rad[0]} and {M.dir[0][2]/rad[1]},
    respectively.  Sets {y[2]} and {y[3]} to the corods of the difference between
    the vector {u = (M.dir[1][1],M.dir[1][2])} and the vector {(1,0)},
    divided by {rad[2]} and {rad[3]}, respectively. */

void hr2_pmap_opt_similarity_decode(double y[], double rad[], hr2_pmap_t *M);
  /* Converts four parameters {y[0..3]} to a similarity map {M}, 
    as the inverse of {hr2_pmap_opt_similarity_encode}.
  
    More precisely, assumes that {y,rad} have at least 4 elements each,
    and {M} is a similarity map with {M.dir[0][0]=M.inv[0][0]=1}. Sets {M.dir[0][1]}
    and {M.dir[0][2]} to {y[0]*rad[0]} and {y[1]*rad[1]}, respectively.
    Then computes the vector {u = (1,0) + (y[2]*rad[2],y[3]*rad[3])},
    stores it into {M.dir[1][1..2]}, and the perpendicular vector 
    {v = (-u[1],+u[0])} into {M.dir[2][1..2]}.  Then
    sets {M.inv} to be the inverse of {M.dir}. */

#endif
