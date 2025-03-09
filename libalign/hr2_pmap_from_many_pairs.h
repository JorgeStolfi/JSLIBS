#ifndef hr2_pmap_from_many_pairs_H
#define hr2_pmap_from_many_pairs_H

/* Tools for optimizing projective maps. */
/* Last edited on 2025-02-16 20:00:50 by stolfi */ 

#include <stdint.h>

#include <bool.h>
#include <r2.h>
#include <hr2.h>

#include <hr2_pmap_opt.h>

hr2_pmap_t hr2_pmap_from_many_pairs
  ( hr2_pmap_type_t type,
    uint32_t np,
    r2_t p1[], 
    r2_t p2[], 
    double w[],
    int32_t maxIter,
    double maxErr,
    bool_t verbose
  );
  /* Computes the projective matrix {M} such that the Cartesian points
    {p1[0..np-1]} mapped by {M} are as close as possible to the
    Cartesian points {p2[0..np-1]}, with importance weights {w[0..np-1]}.  
    
    The points must all be finite and the weights should be all
    non-negative. Effectively, each pair {p1[i],p2[i]} counts as much as
    {w[i]} copies of that pair. If {w} is {NULL}, all points assumed to
    have the same weight.
    
    The /point rank/ of the given point pairs and weights is the size {nr} of
    the largest subset of points pairs {p1[i],p2[i]} that have positive
    weight {w[i]} and such that the points on each side ({p1} or {p2}) are all
    distinct and do not include three collinear points. The optimum map
    is well-defined only if the point rank {nr} is at least
    {nr_min(type) = hr2_pmap_from_many_pairs_required_rank(type)}. Namely, {nr}
    must be at least 0 for identity, 1 for translation, 2 for similarity
    or congruence, 3 for affine, and 4 for general projective.
    
    If the point rank is insufficient for the requested {type}, the
    procedure it will implicilty downgrade the {type} until it becomes
    sufficient. Thus, for example, if {type} requests a general affine
    or projective map, but the actual point rank {nr} is 2, then the
    {type} will be downgraded to {hr2_pmap_type_SIMILARITY}. If the rank
    is only 1, it will be downgraded to translation (skipping
    congruence).
    
    In general, the returned map {M} will take each point {p1[i]} near
    to, but not always onto, the corresponding point {p2[i]}. If the
    requested map {type} is identity, the returned map {M} is in fact
    always the identity, ignoring the points. Otherwise, the procedure
    tries to choose {M} so that it minimizes the weighted average of the
    squared distances as computed by {hr2_pmap_mismatch_sqr(M,np,p1,p2,w)}'.
     
    If the requested map type is a translation, the result is simply the
    translation that moves the barycenter of the points {p1[0..np-1]} to
    the barycenter of the points {p2[0..np-1]}, computed with the
    weights {w[0..np-1]}.
   
    For the other map types, the procedure uses the iterated
    edge-divided simplex method for multivariate optimization, with a
    maximum {maxIter} iterations. The {report} function, if not {NULL},
    is called every time the goal function is evaluated by that method.
    The iteration stops when the procedure thinks that it found a matrix
    that is less than {maxErr} away from optimum. This distance roughly
    refers to the distance between mapped point position. E. g., if the
    point coordinates are measured in units of pixels, {maxErr} should
    typically be a large fraction less than 1.
    
    For the general projective map {type}, an extra term is added to
    the goal function above, to remove the spurious degree of freedom
    due to homogeneous matrix scaling. */

uint32_t hr2_pmap_from_many_pairs_required_rank(hr2_pmap_type_t type);
  /* Returns the minimum number of point pairs that
    {hr2_pmap_from_many_pairs_try_find_frame_pair} must find in order to
    compute the optimum map with the given {type}. Namely, 0 for
    identity, 1 for translation, 2 for similarity or congruence, 3
    for affine, and 4 for general projective. */

#endif
