#ifndef hr2_pmap_opt_congruence_H
#define hr2_pmap_opt_congruence_H

/* Last edited on 2023-10-20 16:21:25 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <hr2.h>

void hr2_pmap_opt_congruence_encode(hr2_pmap_t *M, double rad[], double y[]);
  /* Converts a congruence map {M} to three parameters {y[0..2]}.
  
    More precisely, assumes that {M} is a congruence map (a combined
    rotation and translation) with {M.dir[0][0]=1}, and that {y} and
    {rad} have at least 3 elements each. Sets {y[0]} and {y[1]} to
    {M.dir[0][1]/rad[0]} and {M.dir[0][2]/rad[1]}, respectively. Then
    determines the rotation angle {ang} (in radians) as the angle from
    the {X} axis to the vector {(M.dir[1][1],M.dir[1][2])}, and sets
    {y[2]} to {ang/rad[2]}. */

void hr2_pmap_opt_congruence_decode(double y[], double rad[], hr2_pmap_t *M);
  /* Converts two parameters {y[0..2]} to a congruence map {M}, 
    as the inverse of {hr2_pmap_opt_congruence_encode}.
  
    More precisely, assumes that {y} and {rad} have at least 3 elements each,
    and {M} is a congruence map with {M.dir[0][0]=M.inv[0][0]=1}. Sets {M.dir[0][1]}
    and {M.dir[0][2]} to {y[0]*rad[0]} and {y[1]*rad[1]}, respectively.
    Sets the submatrix {M.dir[1..2][1..2]} to a rotation by angle {ang=y[2]*rad}
    (in radians). Then sets {M.inv} to be the inverse of {M.dir}. */

#endif
