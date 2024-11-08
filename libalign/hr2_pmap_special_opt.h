#ifndef hr2_pmap_special_opt_H
#define hr2_pmap_special_opt_H

/* Tools for optimizing projective maps. */
/* Last edited on 2024-11-08 09:18:21 by stolfi */ 

#define _GNU_SOURCE

#include <bool.h>
#include <stdint.h>
#include <r3x3.h>
#include <i2.h>
#include <hr2.h>
#include <hr2_pmap.h>
#include <hr2_pmap_opt.h>

void hr2_pmap_special_opt_quadratic
  ( hr2_pmap_type_t type,     /* Desired map type. */
    sign_t sgn,               /* Handedness, when applicable. */
    hr2_pmap_opt_func_t *f2,  /* Goal function to minimize. */
    hr2_pmap_opt_pred_t *OK,  /* Client acceptance function. */
    double yrad[],            /* Scaling factor for each encoding element. */
    int32_t maxIter,          /* Max outer loop iterations. */
    hr2_pmap_t *M,            /* (IN/OUT) The projective map to adjust. */
    double *f2M_P,            /* (OUT) Goal function value for the output {*M}. */
    bool_t verbose
  );
  /* Finds a projective map {*M} of the specified {type} and {sgn}
    that minimizes the function {*f2}, using quadratic optimization.   
    The {type} cannot be {hr2_pmap_type_GENERIC} or {hr2_pmap_type_NONE}.
    
    If it is {hr2_pmap_type_IDENTITY}, the procedure simply sets {*M} to
    the identity map or the {XY}-swap map, depending on {sgn}, ignoring
    the other parameters, and does no optimization.
    
    Otherwise, on input, {*M} must contain an initial guess for the
    projective map. On output, it will contain the optimized map. The
    procedure uses up to {maxIter} rounds of the simplex vertex-edge
    (SVE) method, but will stop earlier if {OK} is not {NULL} and it
    finds a map {M} for which {OK(M)} is true. The number of variables
    in the optimization depends on the map {type}.
    
    The {sgn} parameter indicates whether the map should be orientation-preserving
    ({dir = +1}) or orientation-reversing ({dir = -1}).
    
    The parameter {yrad[0..ny-1]} determines maximum change allowed in
    the encoding {y[0..ny-1]} of the map {M}. Its length {ny} should be
    {hr2_pmap_encode_num_parameters(type)}. Ideally the values should be
    such that, for small {d}, a change of {d*yrad[ky]} on the encoding
    element {y[ky]} will have the same effect on {f2(M)}, for any {ky}.
    
    In any case, the output value of {*f2M_P} will be the value of {f2(M)}.
  
    The mismatch function {f2} had better have a single minimum within
    the search region, and preferably be approximately quadratic on all
    variable element of the affine map within that region. */

void hr2_pmap_special_opt_1D_plot
  ( char *outPrefix,
    char *tag,
    char *stage,
    hr2_pmap_t *M,
    hr2_pmap_type_t type,
    sign_t sgn,
    hr2_pmap_opt_func_t *f2,
    double yrad[],
    int32_t nu,
    int32_t ns
  );
   /* Writes to a file called "{outPrefix}-{tag}-{stage}-1D-plot.txt" the values
    {f2(N(ks,ku))}, for {ks} in {-ns..+ns} and {ky} in {0..nu};
    where {N(ks,ku)} is the matrix {M} modified in various ways and amounts.
    
    The modifications depend on the {type}, and happen in the space
    {\RR^ny} of dimension {ny=hr2_pmap_encode_num_parameters(hr2_pmap_type_t
    type)} which is an encoding of all matrices of type {type}.
    
    The procedure chooses {nu} unit directions {u[0..nu-1][0..ny-1]} in that space,
    and lets {N(ks,ku)} be the decoding of {yn(ks,ky)}, where
    {{yn(ks,ku)[ky] = ym[ky] + (ks/nd)*u[ky]*yrad[ky]} for each {ky} in {0..ny-1},
    and {ym[0..ny-1]} is the encoding of {M}.
    
    The matrix {M} should have handedness {sgn}, and the matrices {N(ks,ku)} 
    will have the same {sgn}.
    
    The directions {u} will include up to {nu} cardinal directions of
    {\RR^ny}, and some random diagonal ones if {nu > ny}.
    
    The file will have {2*ns+1} data lines. Each line of the file will
    have the format "{ks} {ds} {f2(N(ks,0))} ... {f2(N(ks,np-1))}" where
    {ks} is a sample index from {-ns} to {+ns} and {ds} is the
    displacement distance {ks*yrad/ns}.
    
    The type must not be {hr2_pmap_type_IDENTITY}, {hr2_pmap_type_GENERIC},
    or {hr2_pmap_type_NONE}. */

#endif
