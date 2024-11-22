#ifndef hr2_pmap_from_many_pairs_aux_H
#define hr2_pmap_from_many_pairs_aux_H

/* Headers of some internal procedures of {hr2_pmap_from_many_pairs.c} for testing. */
/* Last edited on 2024-11-21 21:17:30 by stolfi */ 

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <r2.h>
#include <hr2.h>

#include <hr2_pmap_opt.h>
    
/* INTERNAL PROCEDURES */

void hr2_pmap_from_many_pairs_try_find_frame_pair
  ( hr2_pmap_type_t type_req,
    int32_t np,
    r2_t p1[],
    r2_t p2[],
    double w[],
    int32_t *nr_P,                /* OUT */
    int32_t ixr[]                 /* OUT */
  );
  /* Among the weighted pairs of points 
    {p1[i],p2[i],w[i]}, for {i} in {0..np-1}, tries to find 
    four pairs such that no three of the four points on each side are
    collinear.
    
    If it cannot find such four pairs, tries to find
    three pairs such that the three points on each side are not collinear.
    
    If it still does not succeed, it tries to find at least two pairs
    such that the two points on each side are distinct.
    
    If that too fails, tries to find one pair with positive 
    weight.  
    
    Which too may fail.
    
    In any case, returns in {*nr_P} the number {nr} of points pairs
    found (0 to 4), and in {ixr[0..nr-1]} the indices {i} of those
    pairs.
    
    !!! Currently, if {type} is general projective, the procedure may fail
    to find the 4 required point pairs in some cases when they do exist.
    For instance, when {np=5}, point {p1[3]} is the midpoint of {p1[0]}
    and {p1[1]}, and point {p1[4]} is the midpoint of {p1[0]}
    and {p1[2]}. Fixing this bug does not seem easy, though. !!! */

hr2_pmap_t hr2_pmap_from_many_pairs_initial
  ( hr2_pmap_type_t type,
    int32_t np,
    r2_t p1[],
    r2_t p2[], 
    double w[],
    int32_t nr,
    int32_t ixr[],
    sign_t sgn,
    int32_t class
  );
  /* Returns an initial map {M} of the specified {type} that takes the
    /source frame points/ {p1[ixr[0..nr-1]]} onto or near to 
    the /destination frame points/ {p2[ixr[0..nr-1]]}.  
    
    The number of points {nr} in the should be adequate for the type,
    as specified by {hr2_pmap_from_many_pairs_required_rank}. 
    
    If {type} is the identity, the returned map is the identity, 
    that is obviously optimal.
    
    If {type} is translation, the procedure ignores the parameters {nr,ixr}
    and instead returns the translation that takes the barycenter of {p1[0..np-1]}
    to that of {p2[0..np-1]}; which is optimal too.
    
    If {type} is congruence, the returned map {M} will take the midpoint
    of the points {a1=p1[ixr[0]]} and {b1=p1[ixr[1]]} to the midpoint of
    {a2=p2[ixr[0]]} and {b2=p2[ixr[1]]}. If {type} is similarity, it
    takes {a1} to {a2} and {b1} to {b2}. In either case, the map will be
    such that the direction {M(a1)->M(b1)} will match that of {a2->b2}
    
    If {type} is affine, the returned map {M} will take the 
    three points {a1,b1,c1=p1[ixr[0..2]]} to {a2,b2,c2=p2[ixr[0..2]]},
    or to {a2,c2,b2}. 
    
    If {type} is generic projective, the returned map {M} will be one
    of the maps that take the three hither points with Cartesian
    coordinates {a1,b1,c1} to three points, hither or yonder, with
    Cartesian coordinates {a2,b2,c2} or {a2,c2,b2}.

    The parameter {sgn} is relevant to when type is anything other than
    identity or translation, and specifies that the resulting map should
    be orientation-preserving ({sgn = +1}) or orientation-reversing
    ({sgn = -1}). See the {sgn} parameter of
    {hr2_pmap_congruence_from_point_and_dir} and
    {hr2_pmap_similarity_from_two_points}. For general affine and
    projective maps, setting {sgn = -1} gives a map that takes {b1,c1}
    to {c2,b2} instead of {b2,c2}.
    
    The parameter {class} is relevant for the generic projective map
    type, in which case it should be an integer in {0..3}. The two lower
    bits {class&1} and {class&2} specify whether the hither points with
    Cartesian coordinates {b1} and {c1} are to be mapped to hither or
    yonder points with Cartesian coordinates {b2} and {c2}. See
    {hr2_pmap_r2_from_sign_class}.
    
    If {type} is general projective and {class} is not 0, the resulting map
    will take some hither points to yonder points, some finite points to
    infinite points, and vice-versa.  In all other cases, including
    general projective with {class=0}, the resultimg map will be an affine map
    that takes finite hither points to finite hither points. */

void hr2_pmap_from_many_pairs_optimize
  ( hr2_pmap_type_t type,
    sign_t sgn,
    int32_t np,
    r2_t p1[],
    r2_t p2[],
    double w[],
    int32_t maxIter,
    double maxErr,
    hr2_pmap_t *M,     /* (IN/OUT) The projective map to adjust. */
    double *f2M_P,     /* (IN/OUT) Goal function value for {*A}. */
    bool_t verbose
  );
  /* Uses quadratic optimization to modify the map {M} so as to mimimize
    the goal function {f2(M)} that is the squared distance between the
    points {M.dir(p1[i])} and {p2[i]} and the squared distance between
    {p1[i]} {M.inv(p2[i])}, averaged over all {i} in {0..np-1} with
    weights {w[i]}. (If {w} is {NULL}, assumes that the weights are all
    1.0.)
    
    The optimization assumes that, on entry, {M} is a projective map of
    the specified {type} and {sgn}, which will be preserved, where
    applicable. If {type} is identity or translastion, {sgn} must be 0
    or {+1}, and is ignored; otherwise it must be {-1} or {+1}. On exit,
    {M} will be the optimized map, and {*f2M_P} will be the value of
    {f2(M)}.
    
    The procedure uses the iterated divided-edge simplex method of
    quadratic minimization, with at most {maxIter} iterations, stopping
    when the average discrepancy is {maxErr^2}.  The maximum change
    in any array element is */

#endif
