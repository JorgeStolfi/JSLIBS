#ifndef r2_aff_adjust_H
#define r2_aff_adjust_H

/* Tools for optimizing affine maps. */
/* Last edited on 2020-10-16 03:36:47 by jstolfi */ 

#define _GNU_SOURCE

#include <bool.h>
#include <r2.h>
#include <r2x2.h>
#include <r2_aff_map.h>
#include <i2.h>

typedef double r2_aff_adjust_func_t(r2_aff_map_t *A); 
  /* Type of a procedure that evaluates some badness function of
    the affine map {*A}. The function had better achieve a minimum
    value for some finite {*A}. */

void r2_aff_adjust_quad
  ( r2_aff_adjust_func_t *f2,  /* Goal function to minimize. */
    r2_aff_map_t *R,     /* Max adjustment for {*A}. */
    double tol,          /* Desired relative adjustment precision for {*A}. */
    r2_aff_map_t *A,     /* (IN/OUT) The affine map to adjust. */
    double *f2P          /* (OUT) Goal function value for the output {*A}. */
  );
  /* Similar to {r2_aff_adjust_enum}, but uses quadratic optimization
    instead of enumeration. 

    The meximum adjustment of each element of the affine map {*A} is the
    corresponding element of {*R}. Its accuracy will be about {tol}
    times that maximum adjustment. If the max adjustment is zero, the
    element is assumed to be fixed and is not adjusted.
    
    The mismatch function {f2} had better have a single minimum within
    the search region, and preferably be approximately quadratic on all
    variable element of the affine map within that region. */
 
#endif
