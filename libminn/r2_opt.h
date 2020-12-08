#ifndef r2_opt_H
#define r2_opt_H

/* Tools for optimizing a vector of points on the plane. */
/* Last edited on 2017-06-06 21:40:27 by stolfilocal */ 

#include <bool.h>
#include <r2.h>
#include <i2.h>

typedef double r2_opt_goal_func_t(int ni, r2_t p[], i2_t iscale); 
  /* Type of a function that evaluates the goal function {F} to optimize.
    for points {p[0..ni-1]}, at a specified scale of resolution {iscale}.
    
    The goal function should be computed as if its domain was scaled by
    {1/2^iscale.c[j]} along axis {j}, with each alignment point {p[i]}
    implicitly scaled by that same amount. This scaling is relevant for
    functions that are defined by integrals or averages over some
    neighborhood of the points {p[i]}.
    
    The size of the neighborhood should be the the same for all scales,
    when measured on the scaled function. In other words, the
    neighnorhood should be centered on {p[i]} and should have size
    proportional to {2^iscale.c[j]} along axis {j}, when applied to the
    unscaled function.
    
    For example, the goal function could be some rms discrepancy between
    {ni} images and a reference image, that considers some fixed
    neigborhood of the point {p[i]} in the domain of each image {i}.
    Then, if {iscale.c[j]} is nonzero, each image and the corresponding
    point {p[i]} should be scaled by {1/2^iscale.c[j]} before performing
    the comparison.
    
    The function had better be C2 near the optimum alignment for
    quadratic optimization (see below). */

/* OPTIMIZATION FUNCTIONS

  Each procedure in this section adjusts a point vector {p[0..ni-1]} so
  as to minimize the goal function function {f2(ni,p,iscale)}.

  On input, {p[0..ni-1]} must be a guess for the optimum point vector.
  On output, {p[0..ni-1]} will be an improved point vector, which
  differs from the given one by an integer multiple of {astp[i].c[j]} in
  each axis {j=0..1}, not exceeding {arad[i].c[j]} in absolute value.

  If {arad[i].c[j]} is positive, {astp[i].c[j]} must be positive. If
  {arad[i].c[j]} is zero or less than {astp[i].c[j]}, the coordinate
  {p[i].c[j]} is considered fixed, and not changed.

  The value of {f2(ni,p,iscale)} is returned on {*f2p}. */

    
void r2_opt_single_scale_enum
  ( int ni,                   /* Number of points to optimize. */
    i2_t iscale,              /* Object scaling exponent along each axis. */  
    r2_opt_goal_func_t *f2,   /* Function that evaluates the goal function. */
    r2_t arad[],              /* Max coordinate adjustment for each point. */
    r2_t astp[],              /* Adjustment step for each point. */
    r2_t p[],                 /* (IN/OUT) Corresponding points in each object. */
    double *f2p,              /* (OUT) Goal function value at the best solution. */
    bool_t debug              /* Prints diagnostics if true. */
  );
  /* Adjusts a point vector {p[0..ni-1]} so as to minimize the goal
    function function {f2(ni,p,iscale)}. Uses exhaustive enumeration of
    all valid point vectors. */

void r2_opt_single_scale_quadopt
  ( int ni,                   /* Number of points to optimize. */
    i2_t iscale,              /* Object scaling exponent along each axis. */  
    r2_opt_goal_func_t *f2,   /* Function that evaluates the goal function. */
    r2_t arad[],              /* Max coordinate adjustment for each point. */
    r2_t astp[],              /* Desired adjustment precision for each point. */
    r2_t p[],                 /* (IN) Initial guesses; (OUT) Best solution. */
    double *f2p,              /* (OUT) Goal function value at the best solution. */
    bool_t debug              /* Prints diagnostics if true. */
  );
  /* Adjusts a point vector {p[0..ni-1]} so as to minimize the goal
    function function {f2(ni,p,iscale)}. Uses iterated quadratic
    minimization.
    
    The goal function function {f2} had better have a single minimum within
    the search region, and preferably be approximately quadratic on {p}
    within that region. */

void r2_opt_multi_scale
  ( int ni,                  /* Number of points to optimize. */
    r2_opt_goal_func_t *f2,  /* Function that evaluates the goal function. */
    bool_t quadopt,          /* Use quadratic optimization? */
    r2_t arad[],             /* Max coordinate adjustment for each point. */
    r2_t astp[],             /* Adjustment step or desired precision for each point. */
    r2_t p[],                /* (IN) Initial guesses; (OUT) Best solution. */
    double *f2p,              /* (OUT) Goal function value at the best solution. */
    bool_t debug              /* Prints diagnostics if true. */
  );
  /* Adjusts a point vector {p[0..ni-1]} so as to minimize the goal
    function function {f2(ni,p,(0,0))}. Uses a multiscale point vector
    strategy, with variable step sizes, to efficiently search large
    coordinates {p[i].c[j]} with large uncertainty (large ratio
    {arad[i].c[j]/astp[i].c[j]}).
    
    The procedure performs an enumerative or quadratic optimization of
    the adjustment vector at various scales {iscale}, starting with
    some coarsest scale {mscale} and ending with scale {(0,0)}.
    Successive scales have at least one of the coordinates reduced by 1.
    
    At each scale, the procedure uses {r2_opt_single_scale_enum}
    or {r2_opt_single_scale_quadopt}, 
    with modified search steps {stp} and modified search radii {rad}
    such that (1) {stp[i].c[j]} is {astp[i].c[j]*2^iscale.c[j]}, or
    zero; and (2) {rad[i].c[j]} is about {2*stp[i].c[j]}, or zero. Thus,
    at each scale the adjustment range for {p[i].c[j]} is only a small
    multiple of the adjustment step/precision
    
    The procedure starts the search at a scale {mscale} such that
    {mscale.c[j]} is the smallest integer such that
    {astp[i].c[j]*2^iscale.c[j] > arad[i].c[j]}, for all {i} such that
    {arad[i].c[j] > 0}. */

double r2_opt_rel_disp_sqr(int ni, r2_t p[], r2_t q[], r2_t arad[], r2_t astp[]);
  /* Computes the total squared displacement between {p[0..ni-1]} and 
    {q[0..ni-1]} relative to the radius {arad[i]}, that is,
    
      { SUM { (p[i].c[j] - q[i].c[j])/arad[i].c[j]|^2 : i=0..ni-1,j=0..1 } }
    
    except that terms where {arad[i].c[j]} is zero or less than
    {astp[i].c[j]} are ignored. Clients may want to add an appropriate
    multiple of this function to the goal function in order to bias the
    search towards the neighborhood of the initial guess. */

#endif
