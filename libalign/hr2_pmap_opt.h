#ifndef hr2_pmap_opt_H
#define hr2_pmap_opt_H

/* Tools for optimizing projective maps. */
/* Last edited on 2023-10-08 11:37:24 by stolfi */ 

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

hr2_pmap_t hr2_pmap_from_many_pairs
  ( int32_t np,
    r2_t p1[], 
    r2_t h1[],
    r2_t p2[], 
    r2_t h2[],
    hr2_pmat_type_t type,
    int32_t maxIter,
    double maxErr,
    hr2_pmap_opt_report_proc_t report
  );
  /* Computes the projective matrix {M} such that the Cartesian points
    {p1[0..np-1]} mapped TWICE by {M} are as close as possible to the
    Cartesian points {p2[0..np-1]}. The number of points should preferably be 
    the minimum needed by the type, or more:
    
    (If {np} is 4, one should probably use {hr2_pmap_from_four_points}
    instead.)
    
    Uses the iterated edge-divided simplex method for optimization, with
    a maximum {maxIter} iterations. Stops when it thinks that it found a
    matrix that is less than {maxErr} away from optimum. This distance
    roughly refers to the distance between mapped point position. E. g.,
    if the points are pixel coordinates in two images, {maxErr} should
    typically be a large fraction less than 1.
    
    More precisely, the procedure tries to minimize the squared distance
    between each {p1[i]} mapped by {M} and {p2[i]} mapped by {M^{-1}},
    averaged over all {i}. Thus {M} and {M^{-1}} bring the points to a
    common ``middle ground'' plane where they are most similar.
    
    In general there are many maps that would do that, because the
    square root of a matrix is not unique. Also, the map {M^2} is not
    unique if each set of points lies on a common quadric curve. The
    arguments {h1} and {h2} are meant to solve that problem.
    
    The arrays {h1} and {h2} should have four elements each, which are
    taken to be the corners of quadrilaterals {Q1} and {Q2} on the
    domain of point lists {p1} and {p2}, respectively. Typically, but
    not necessarily, each quadrilateral will enclose the corresponding
    point list. The procedure is lightly biased towards the map {M} such
    that {M.dir} causes the smallest deformation of {Q1}, and {M.inv}
    causes the smallest deformation of {Q2}. If {h1} is {NULL}, the
    bounding rectangle of the points {p1} is used instead; and similarly
    for {h2} and {p2}.
    
    The {report} function, if not {NULL}, is called every time the 
    goal function is evaluated by that method. */

void hr2_pmap_opt_quad
  ( hr2_pmap_opt_func_t *f2,  /* Goal function to minimize. */
    r3x3_t *R,                   /* Max adjustment for {A.dir}. */
    double rtol,                 /* Desired relative adjustment precision for {*A}. */
    hr2_pmap_t *A,               /* (IN/OUT) The affine map to adjust. */
    double *f2AP                 /* (OUT) Goal function value for the output {*A}. */
  );
  /* Finds a projective map {*A} that minimizes the function {*f2},
    using quadratic optimization. 
    
    On input, {*A} must have the initial guess for the projective map.
    On output, it will have the optimized map, and {*f2AP} will be
    the value of {f2(A)}.

    The meximum adjustment of each element of the direct map {A.dir} is the
    corresponding element of {*R}. Its accuracy will be about {tol}
    times that maximum adjustment. If the max adjustment is zero, the
    element is assumed to be fixed and is not adjusted.
    
    Because of the homogeneous scaling invariance, is an error if all 9
    elements of {A.dir} are variable (that is, all 9 elements of {*R} are
    nonzero). At least one of them must be fixed.
    
    The mismatch function {f2} had better have a single minimum within
    the search region, and preferably be approximately quadratic on all
    variable element of the affine map within that region. */
 
void hr2_pmap_opt_get_var_elems
  ( r3x3_t *M, 
    r3x3_t *R, 
    double *Mp[],
    double Re[],
    int32_t *nvP
  );
  /* Identifies the elements of {*M} that are to be adjusted 
    (that is, whose corresponding element of {*R} is nonzero.
    Returns the number {nv} of such elements (a number in {0..8}) in {*nvP},
    the addresses of those elemens of {*M} in {Mp[0..nv-1]},
    and the values of the corresponding elements of {*R}
    in {Re[0..nv-1]}.  These vectors should have at leat 8 elements.
    The procedure fails if all 9 elements of {*R} are nonzero.*/

#endif
