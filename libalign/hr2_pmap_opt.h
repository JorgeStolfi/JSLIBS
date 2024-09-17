#ifndef hr2_pmap_opt_H
#define hr2_pmap_opt_H

/* Tools for optimizing projective maps. */
/* Last edited on 2024-09-17 19:58:14 by stolfi */ 

#define _GNU_SOURCE

#include <bool.h>
#include <stdint.h>
#include <r3x3.h>
#include <hr2.h>
#include <hr2_pmap.h>
#include <i2.h>

typedef double hr2_pmap_opt_func_t(hr2_pmap_t *M); 
  /* Type of a procedure that evaluates some badness function of the
    projective map {*M}. The function had better achieve a minimum value
    for some finite {*M}. */

void hr2_pmap_opt_quadratic
  ( hr2_pmap_type_t type,     /* Desired map type. */
    sign_t sgn,               /* Handedness, when applicable. */
    hr2_pmap_opt_func_t *f2,  /* Goal function to minimize. */
    int32_t maxIter,          /* Max outer loop iterations. */
    double maxMod,            /* Maximum change allowed in each matrix element. */
    double f2Stop,            /* Stopping criterion. */
    hr2_pmap_t *M,            /* (IN/OUT) The projective map to adjust. */
    double *f2M_P,            /* (OUT) Goal function value for the output {*M}. */
    bool_t verbose
  );
  /* Finds a projective map {*M} of the specified {type}
    that minimizes the function {*f2}, using quadratic optimization.   
    The {type} cannot be {hr2_pmap_type_NONE}.
    
    If it is {hr2_pmap_type_IDENTITY}, the procedure simply sets {*M} to
    the identity map, ignoring the other parameters, and does no
    optimization.
    
    Otherwise, on input, {*M} must contain an initial guess for the
    projective map. On output, it will contain the optimized map. The
    procedure uses up to {maxIter} rounds of the simplex vertex-edge
    (SVE) method, but will stop earlier if it finds a map {M} for which
    {f2(M) <= f2Stop}. The number of variables in the optimization
    depends on the map {type}.
    
    The {sgn} parameter indicates whether the map should be orientation-preserving
    ({dir = +1}) or orientation-reversing ({dir = -1}).  It is ignored if the {type}
    is identity or translation, which are always orientation-preserving.
    
    The parameter {maxMod} determines maximum change allowed in the
    matrix {M.dir}.
    
    In any case, the output value of {*f2M_P} will be the value of {f2(M)}.
  
    The mismatch function {f2} had better have a single minimum within
    the search region, and preferably be approximately quadratic on all
    variable element of the affine map within that region. */

void hr2_pmap_opt_1D_plot
  ( char *outPrefix,
    char *tag,
    hr2_pmap_t *M,
    hr2_pmap_type_t type,
    sign_t sgn,
    hr2_pmap_opt_func_t *f2,
    int32_t nu,
    double urad,
    int32_t ns
  );
   /* Writes a file called "{outPrefix}-{tag}-1D-plot.txt" the values
    {f2(N(ks,ku))}, for {ks} in {-ns..+ns} and {ky} in {0..nu};
    where {N(ks,ku)} is the matrix {M} modified in various ways and amounts.
    
    The modifications depend on the {type}, and happen in the space
    {\RR^nz} of dimension {nz=hr2_pmap_encode_num_parameters(hr2_pmap_type_t
    type)} which is an encoding of all matrices of type {type}.
    
    The procedure chooses {nu} directions {u[0..nu-1]} in that space,
    and lets {N(ks,ku)} be {M} displaced by {(ks*urad/ns)*u[ku]} in that space.
    
    If {type} is not translation, the matrix {M} should have handedness
    {sgn}, and the matrices {N(ks,ku)} will have the same {sgn}.
    
    The directions {u} will include up to {nu} cardinal directions of
    {\RR^nz}, and some random diagonal ones if {nu > nz}.
    
    The file will have {2*ns+1} data lines. Each line of the file will
    have the format "{ks} {ds} {f2(N(ks,0))} ... {f2(N(ks,np-1))}" where
    {ks} is a sample index from {-ns} to {+ns} and {ds} is the
    displacement distance {ks*urad/ns}.
    
    The type must not be {hr2_pmap_type_IDENTITY}. If it is
    {hr2_pmap_type_GENERIC}, the plotted value includes a small
    multiple of bias term that is designed to remove the spurious degee
    of freedom in {\RR^nz} that consists of homogeneous scaling of the
    whole map matrices. */

#endif
