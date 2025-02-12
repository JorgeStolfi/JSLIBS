#ifndef r2_align_quadopt_H
#define r2_align_QUadopt_H

/* Tools for optimizing translational alignment of 2D objects by iterated quadratic optimization. */
/* Last edited on 2024-12-05 10:20:14 by stolfi */ 

#include <stdint.h>

#include <bool.h>
#include <r2.h>
#include <i2.h>
#include <r2_align.h>

void r2_align_quadopt
  ( int32_t ni,               /* Number of objects to align. */
    r2_align_mismatch_t *F2,  /* Function that evaluates the mismatch between the objects. */
    r2_t arad[],              /* Max delta vector coordinates for each object. */
    bool_t bal,               /* True if alignment vector adjustments should be balanced. */
    double tol,               /* Desired precision. */
    r2_t p[],                 /* (IN/OUT) Initial and optimized alignment vector. */
    double *F2val_P           /* (OUT) Mismatch for the computed alignment vector. */
  );
  /* Adjusts an alignment vector {p[0..ni-1]} for certain {ni} objects so
    as to minimize a given mismatch function {F2}.
    
    On input, {p[0..ni-1]}, must be a guess for the optimum alignment.
    On output, {p[0..ni-1]} will be the best alignment found. The value
    of {F2(ni,p)} is returned on {*F2val_P}.
    
    Uses iterated quadratic minimization within the search ellipsoid
    {\RD} centered at the initial alignment vetor {p} and defined by the
    parameters {arad,bal}. The search stops when the procedure believes
    that it is within distance {tol} from the optimum.
    
    The mismatch function {f2} had better have a single minimum within
    the search region, and preferably be approximately quadratic on {p}
    within that region. */

#endif
