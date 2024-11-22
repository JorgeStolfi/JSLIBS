/* Affine maps in 2D oriented projective plane {\RT^2}. */
/* Last edited on 2024-11-21 20:20:27 by stolfi */ 
   
#ifndef hr2_pmap_affine_H
#define hr2_pmap_affine_H

#define _GNU_SOURCE
#include <stdint.h>

#include <sign.h>
#include <bool.h>
#include <r2.h>
#include <r2x2.h>
#include <hr2_pmap.h>

double hr2_pmap_affine_discr_sqr(hr2_pmap_t *M, hr2_pmap_t *N);
  /* Assumes that {M} and {N} are affine maps.  Returns the total
    squared mismatch between them, defined as 
    
      { INTEGRAL { |(A - B)(u)|^2 du: u \in \RS^1 } }
    
    */

hr2_pmap_t hr2_pmap_affine_from_mat_and_disp(r2x2_t *E, r2_t *d);
  /* Returns an affine map {M} that performs the linear map 
    described by the matrix {E} followed by translation by {d}.
    That is, maps the point with Cartesian coordinates {p \in \RR^2} 
    to {p*E + d}. The map will have unit weight (that is, 
    {M.dir.c[0][0] == 1}). */
     
hr2_pmap_t hr2_pmap_affine_from_three_points(r2_t *o, r2_t *p, r2_t *q);
  /* Returns a projective map that takes the points {(0,0)},
    {(0,1)}, and {(1,0)} to {o}, {p}, and {q}, respectively.
    
    The procedure returns a degenerate (non-invertble) map
    if {o,p,q} are collinear.  */

#endif
