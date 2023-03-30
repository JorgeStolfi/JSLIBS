#ifndef r2_align_multiscale_H
#define r2_align_multiscale_H

/* Tools for optimizing a vector of points on the plane. */
/* Last edited on 2023-03-22 19:48:56 by stolfi */ 

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <r2.h>
#include <i2.h>

/* 
  
  ALIGNMENT VECTOR
  
  The procedures in this interface look the two-dimensional alignment of
  {ni} objects. See {r2_align.h} for the concepts of /alignment vector/.
  As in that interface, we denote by {\RC = (\RR^2)^ni} the set of all
  alignment vectors, and by {\BX} and {\BY} the alignment vectors such
  that {\BX[i]} = (1,0)} sand {\BY[i] = (0,1)} for all {i} in {0..ni-1}.

  See that interface also for the definition of the /basic ellipsoid/
  {\RE \sub \RC}, defined by alignment vectors {ctr[0..ni-1]} (the
  /domain center/) and {arad[0..ni-1]} (the /radius vector/); of its
  affine span {\RB}; and of the space of /conformal delta vectors/
  {\RV}, namely the linear subspace {\RV} of {\RC}
  that is tangent to {\RB}.
  
  Note the dimension of {\RB} and {\RV} is {nv = 2*ni - nz}, where
  {nz} is the number of coordinates of {arad} that are zero.
  
  See that interface also for the definition of /balanced conformal delta vector/, and the
  space {\RU \subseteq \RV} of all such delta vectors; of /balanced alignment/,
  and the affine space {\CB = ctr + \RU}; and of the balanced the /search ellipsoid/
  {\CF = \CE \cap \CB}
  
  SEARCH PRECISION
  
  The search is always terminated when the procedure believes that it
  has found an alignment that is within distance {tol} from the] optimum {q}
  in the search region. Namely, such that the sum of the squares 
  of {(p[i].c[j] - q[i].c[j])/tol} is at most 1.

  BALANCED SEARCH
  
  The procedures will also limit the search to alignment vectors 
  that are /balances/ namely such that the differences {p[i].c[j] - p0[i].c[j]}
  add to zero. 
  
  Therefore, the search must have at least two coordinates that are not 
  fixed. Otherwise the search will return the initial guess {p0}
  with no adjustment. 
  */

typedef double r2_align_multiscale_mismatch_t(int32_t ni, r2_t p[], i2_t iscale); 
  /* Type of a function that evaluates the mismatch of {ni} objects.
    in some neighborhood of points {p[0..ni-1]}, at a specified
    scale of resolution {iscale}.  The size of the neighborhood
    is client-defined, but is assumed to be proportional
    to {2^iscale.c[j]} along each axis {j}.
    
    Nevertheless, the function should ideally take about the same time
    whatever the {iscale}. That is, each object should be sampled with a
    fixed number of sampling points, spread over an interval
    proportional to {2^iscale.c[j]} along each axis {j}.
    
    The size of the neighborhood, measured on the scaled object, should
    be the same. That is, measured on the unscaled object, the
    neighborhod should be centered on {p[i]} and have size.
    
    The client must ensure that the sampling is reliable, e.g. by
    interpolating, smoothing, or shrinking the object according to
    Nyquist's criterion.
    
    Said another way, the mismatch should be computed as if every object
    is expanded or reduced by the scale factor {1/2^iscale.c[j]} along
    axis {j}, with a suitable smoothing filter, with each alignment
    point coordinate {p[i].c[j]} implicitly scaled by that same amount;
    and then sampled at a fixed number of points over a neighborhood of
    fixed size, independent of {iscale}.
    
    The function had better be C2 near the optimum alignment for
    quadratic optimization (see below). */

void r2_align_multiscale_single_scale_enum
  ( int32_t ni,                         /* Number of objects to align. */
    i2_t iscale,                        /* Object scaling exponent along each axis. */  
    r2_align_multiscale_mismatch_t *F2, /* Function that evaluates the mismatch between the objects. */
    r2_t arad[],                        /* Max delta vector coordinates for each object. */
    double tol,                         /* Desired precision. */
    r2_t p[],                           /* (IN/OUT) Corresponding points in each object. */
    double *F2val_P                     /* (OUT) Mismatch for the computed alignment vector. */
  );
  /* Adjusts an alignment vector {p[0..ni-1]} for some {ni} objects so
    that they are as similar as possible in the neighborhoords of those
    points, as defined by the mismatch function {F2(ni,p,iscale)}.
    Uses exhaustive enumeration of balanced alignment vectors withing the 
    allowed region defined by {arad} and {tol}.
    
    On input, {p[0..ni-1]}, must be a guess for the optimum alignment.
    On output, {p[0..ni-1]} will be the best alignment found.
    
    The value of {F2(ni,iscale,p)} is returned on {*F2val_P}. */

void r2_align_multiscale_single_scale_quadopt
  ( int32_t ni,                         /* Number of objects to align. */
    i2_t iscale,                        /* Object scaling exponent along each axis. */  
    r2_align_multiscale_mismatch_t *F2, /* Function that evaluates the mismatch between the objects. */
    r2_t arad[],                        /* Max delta vector coordinates for each object. */
    double tol,                         /* Desired precision. */
    r2_t p[],                           /* (IN/OUT) Corresponding points in each object. */
    double *F2val_P                     /* (OUT) Mismatch for the computed alignment vector. */
  );
  /* Adjusts an alignment vector {p[0..ni-1]} for certain {ni} objects so
    that they are as similar as possible in the neighborhoords of those
    points, as defined by the mismatch function {F2(ni,iscale,p)}.
    Uses iterated quadratic minimization within the set of balanced alignment
    vectors allowed by {arad} and {tol}. 
    
    On input, {p[0..ni-1]} must be a guess for the optimum alignment
    vector. On output, {p[0..ni-1]} will be the best alignment found.
    
    The value of {F2(ni,iscale,p)} is returned on {*F2val_P}.
    
    The mismatch function {F2} had better have a single minimum within
    the search region, and preferably be approximately quadratic on {p}
    within that region. */

void r2_align_multiscale
  ( int32_t ni,                         /* Number of objects to align. */
    r2_align_multiscale_mismatch_t *F2, /* Function that evaluates the mismatch between the objects. */
    bool_t quadopt,                     /* Use quadratic optimization? */
    r2_t arad[],                        /* Max delta vector coordinates for each object. */
    double tol,                         /* Desired precision. */
    r2_t p[],                           /* (IN/OUT) Corresponding points in each object. */
    double *F2val_P                     /* (OUT) Mismatch for the computed alignment vector. */
  );
  /* Adjusts an alignment vector {p[0..ni-1]} for certain {ni} objects
    so that they are as similar as possible in the neighborhoords of
    those points, as defined by the mismatch function {F2(ni,(0,0),p)}. Uses
    a multiscale alignment strategy to efficiently search large alignment ranges.
    
    On input, {p[0..ni-1]} must be a guess for the optimum alignment
    vector at scale {(0,0)}. On output, {p[0..ni-1]} will be the best 
    alignment found at scale {(0,0)}.
    
    The value of {F2(ni,(0,0),p)} is returned on {*F2val_P}.
    
    The procedure performs an enumerative or quadratic optimization of
    the delta vector vector at various scales {iscale}, starting with
    some coarsest scale {ismax} and ending with scale {(0,0)}.
    Successive scales have at least one of the coordinates reduced by 1.
    The radii and tolerances are adjusted so that the procedure does 
    aproximately constant work at each scale.
    
    At each scale, the procedure uses
    {r2_align_multiscale_single_scale_enum} or
    {r2_align_multiscale_single_scale_quadopt}, accordng to the parameter
    {quadopt}. */

double r2_align_multiscale_rel_disp_sqr(int32_t ni, r2_t p[], r2_t q[], r2_t arad[]);
  /* Computes the total squared displacement between {p[0..ni-1]} and 
    {q[0..ni-1]} relative to the radius {arad[i]}, that is,
    
      { SUM { (p[i].c[j] - q[i].c[j])/arad[i].c[j]|^2 : i=0..ni-1,j=0..1 } }
    
    except that terms where {arad[i].c[j]} is zero are ignored. Clients
    may want to add an appropriate multple of this function to the goal
    function in order to bias the search towards the neighborhood of the
    initial guess. */

#endif
