#ifndef hr2_pmap_generic_opt_H
#define hr2_pmap_generic_opt_H

/* Tools for optimizing projective maps. */
/* Last edited on 2024-12-05 10:19:40 by stolfi */ 

#include <bool.h>
#include <stdint.h>
#include <r3x3.h>
#include <i2.h>
#include <hr2.h>
#include <hr2_pmap.h>
#include <hr2_pmap_opt.h>

void hr2_pmap_generic_opt_quadratic
  ( sign_t sgn,               /* Handedness, when applicable. */
    hr2_pmap_opt_func_t *f2,  /* Goal function to minimize. */
    hr2_pmap_opt_pred_t *OK,  /* Client acceptance function. */
    int32_t maxIter,          /* Max outer loop iterations. */
    double maxMod,            /* Maximum change allowed in each matrix element. */
    hr2_pmap_t *M,            /* (IN/OUT) The projective map to adjust. */
    double *f2M_P,            /* (OUT) Goal function value for the output {*M}. */
    bool_t verbose
  );
  /* Finds a generic projective map {*M} of the specified sign 
    that minimizes the function {*f2}, using quadratic optimization. 
    
    On input, {*M} must contain an initial guess for the
    projective map. On output, it will contain the optimized map. The
    procedure uses up to {maxIter} rounds of the simplex vertex-edge
    (SVE) method, but will stop earlier if {OK} is not {NULL} and it
    finds a map {M} for which {OK(M)} is true. The number of variables
    in the optimization depends on the map {type}.
    
    The {sgn} parameter indicates whether the map should be orientation-preserving
    ({dir = +1}) or orientation-reversing ({dir = -1}).
    
    The parameter {maxMod} determines maximum change allowed in the
    matrix {M.dir}.
    
    In any case, the output value of {*f2M_P} will be the value of {f2(M)}.
  
    The mismatch function {f2} had better have a single minimum within
    the search region, and preferably be approximately quadratic on all
    variable element of the affine map within that region. */

void hr2_pmap_generic_opt_1D_plot
  ( char *outPrefix,
    char *tag,
    char *stage,
    hr2_pmap_t *M,
    sign_t sgn,
    hr2_pmap_opt_func_t *f2,
    int32_t nu,
    double urad,
    int32_t ns
  );
   /* Writes to a file called "{outPrefix}-{tag}-{stage}-1D-plot.txt" the values
    {f2(N(ks,ku))}, for {ks} in {-ns..+ns} and {ku} in {0..nu}; where
    {nu <= 7} and {N(ks,ku)} is the matrix {M} modified in various ways and
    amounts.
    
    The modifications happen in the space {\RR^{3\times 3}} of real 3x3
    matrices. For each {ku} in {0..nu-1}, the matrix {N(0,ku)} is
    {M.dir} normalized to unit root-mean-square norm, so that it becomes
    a point of the sphere {\RS^8}. The procedure chooses a unit
    direction {u[ku]} in {\RR^{3\times 3}}] that is tangent to {\RS^8}
    at that point, and lets {N(ks,ku)} be {NS(0,ku)} displaced by
    {(ks*urad/ns)*u[ku]} in that space. If {N(ks,ku)} has zero
    determinant or the wrong sign, that line is omitted from the file.
    
    The file will have up to {2*ns+1} data lines. Each line of the file will
    have the format "{ks} {ds} {f2(N(ks,0))} ... {f2(N(ks,np-1))}" where
    {ks} is a sample index from {-ns} to {+ns} and {ds} is the
    displacement distance {ks*urad/ns}. */

#endif
