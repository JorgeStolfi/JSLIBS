#ifndef hr2_pmap_opt_H
#define hr2_pmap_opt_H

/* Tools for optimizing projective maps. */
/* Last edited on 2023-10-20 17:53:57 by stolfi */ 

#define _GNU_SOURCE

#include <bool.h>
#include <stdint.h>
#include <r3x3.h>
#include <hr2.h>
#include <i2.h>

typedef double hr2_pmap_opt_func_t(hr2_pmap_t *A); 
  /* Type of a procedure that evaluates some badness function of the
    projective map {*A}. The function had better achieve a minimum value
    for some finite {*A}. */
    
typedef void hr2_pmap_opt_report_proc_t (hr2_pmap_t *M, double F);
  /* Type of a procedure used by projective map
    optimization functions to report
    the probes made during optimization.  
    
    It is called after every call to the internal goal function. The map {M}
    is the projective map that was tried, and {F} is the value of
    the goal function for it. */

void hr2_pmap_opt_aff_quadratic
  ( hr2_pmap_opt_func_t *f2,  /* Goal function to minimize. */
    hr2_pmat_type_t type,     /* Desired map type. */
    int32_t max_iter,         /* Max outer loop iterations. */
    double max_f2,            /* Stopping criterion. */
    double max_mod,           /* Maximum change allowed in each matrix element. */
    hr2_pmap_t *A,            /* (IN/OUT) The affine map to adjust. */
    double *f2A_P             /* (OUT) Goal function value for the output {*A}. */
  );
  /* Finds a projective map {*A} of the specified {type}
    that minimizes the function {*f2}, using quadratic optimization.   
    The {type} cannot be {hr2_pmat_type_PROJECTIVE};
    it can be {hr2_pmat_type_AFFINE} or any sub-class
    thereof.
    
    On input, {*A} must have the initial guess for the projective map.
    On output, it will have the optimized map, and {*f2A_P} will be
    the value of {f2(A)}.
    
    The procedure uses up to {max_iter} rounds of the simplex vertex-edge
    (SVE) method, but will stop earlier if it finds a map {A} for 
    which {f2(A) <= max_f2}. The number of variables in the optimization 
    depends on the map {type}. 
    
    The parameter {max_mod} determines maximum change allowed
    in the matrix {A}.  The precise interpretation depends on the 
    map type, but generally {A} will the the initial map {A0}
    multiplied by a map {M} of the same {type} such that
    the elements of {M.dir-I} and {M.inv-I} are at most 
    {max_mod} in absolute value, or thereabouts.
  
    The mismatch function {f2} had better have a single minimum within
    the search region, and preferably be approximately quadratic on all
    variable element of the affine map within that region. */

#endif
