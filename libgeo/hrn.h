/* Oriented projective geometry in {N} dimensions.  */
/* Last edited on 2024-11-20 12:10:19 by stolfi */

#ifndef hrn_H
#define hrn_H

#define _GNU_SOURCE
#include <stdint.h>

#include <rn.h>
#include <rmxn.h>
#include <sign.h>

/* Based on HR2.i3, created 1994-05-04 by J. Stolfi. */

/* ORIENTED PROJECTIVE GEOMETRY 

  The /{n}-dimensional oriented projective space/ {T^n} is described in
  
    Jorge Stolfi
    Oriented Projective Geometry: A Framework for Geometric Computations
    237 pp. Academic Press, 1991. ISBN 0-12-672025-8.
  
  POINTS
  
  A point {x} of {T^n} is represented by a vector of {n+1}
  double-precision /homogeneous coordinates/ {[x[0],..x[n]]}, not all
  zero.
  
  Two such vectors represent the same point if and only if one is a
  positive multiple of the other. The point
  {[-x[0],..-x[n]]} is the /antipode/ of {x}, denoted by {¬x}.
 
  HYPERPLANES
  
  A /hyperplane/ is an {(n-1)}-dimensional projective subspace of {T^n}.
  A hyperplane {h} is represented by a vector of {n+1}
  double-precision /homogeneous coeeficients/ <h[0],..h[n]>, not all zero.
  
  Two such vectors represent the same hyperplane if and only if 
  one is a positive multiple of the other. The hyperplane
  {<-h[0],..-h[n]>} is the /opposite/ of {h}, denoted by {¬h}.
  
  DUALITY
  
  The /standard duality/ of {T^n} maps a point {[c[0],..c[n]} to
  the hyperplane {<c[0],..c[n]>}, and vice-versa.
  
  INVALID POINTS AND PLANES

  The homogeneous coordinate vector {[0,0,...,0]} is not a point of
  {T^n}, but is returned sometimes by procedures that normally returns
  a point. It will be called the /invalid point/ and denoted by
  {[0*]}.
  
  Dually, the homogeneous coefficient vector {<0,0,...,0>} is not a
  hyperplane of {T^n}, but is returned sometimes by procedures that
  normally returns a hyperplane. It will be called the /invalid
  hyperplane/ and denoted by {<0*>}.
  
  CANONICAL CARTESIAN INCLUSION
  
  The /canonical inclusion/ of {R^n} in {T^n} is the map that takes
  the point of {R^n} with Cartesian coordinates {(X[0],..X[n-1])} to
  the point of {T^n} with homogeneous coordinates {[x[0],..x[n]] =
  [1,X[0],..X[n-1]]}  (or, more generally, to {[w,w*X[0],.. w*X[n-1]]},
  for any poisitive {w}).

  This map takes affine hyperplanes of {R^n} to hyperplanes of {T^n}.
  It preserves the Euclidean distance, as measured by {rn_dist} in
  {R^n} and {hrn_pt_pt_dist} in {T^n}.
  
  Let's define the /weight/ of a point {x} of {T^n} as its first
  homogeneous coordinate {x[0]}. The range of the canonical inclusion
  map (that is, the canonical isometric model of {R^n} in {T^n}) is
  the /hither part/ of {T^n}, the set of all points with positive
  weight.
  
  The points of {T^n} with negative weight constitute its /yonder part/,
  which is a second isometric model of {R^n} inside {T^n}.
  
  The hither and yonder parts are separated by the /hyperplane at
  infinity/, with coefficients {<1,0,0,...,0>}, consisting of all
  points of {T^n} with zero weight. */
  
void rn_to_hrn(uint32_t n, double P[], double w, double p[]);
  /* Given the Cartesian coordinates {P[0..n-1]} of a point {P} of
    {R^n}, returns the homogeneous coordinates {p[0..n]} of the
    corresponding point {p} of {T^n} with weight {w}, namely {[w,w*P[0],..
    w*P[n-1]]}.
    
    If {w} is positive, the resulting point {p} will be the canonical
    inclusion of point {P} in the ``hither'' part of {T^n}. If {w} is
    negative, the result is the antipode of the canonical inclusion of {P}.
    If {w} is zero, the result is the invalid point {[0*]}. */
    
void hrn_to_rn(uint32_t n, double p[], double P[]);
  /* Given the homogeneous coordinates {p[0..n]} of a point {p}
    of {T^n}, returns the Cartesian coordinates {P[0..n-1]} of the 
    corresponding point {P} of {R^n}, namely {(p[1]/p[0],p[2]/p[0],...,p[n]/p[0])}.
    
    If {p} is in the hither part of {T^n} (that is, {p[0] > 0}), the
    point {P} is the image of {p} under the inverse of the canonical
    inclusion map. If {p} is in the yonder part ({p[0] < 0}), then {P}
    is the inverse image of the antipode of {p}. If {p} is at infinity
    ({p[0] = 0}), the result {P} is undefined (infinite or NaN). */

sign_t hrn_side(uint32_t n, double p[], double h[]);
  /* Returns sign of point {p} relative to hyperplane {h}: 0 on the
    hyperplane, +1 in positive halfspace, -1 in negative halfspace.
    May give inconsistent results for points very close to the
    hyperplane. */

/* PROJECTIVE MAPS */

typedef struct hrn_pmap_t { uint32_t m; uint32_t n; double *dir; double *inv; } hrn_pmap_t;
  /* A projective function from {T^m} to {T^n}.
  
  Field {dir} is the function's matrix, with {m+1} rows and {n+1} columns.
  Field {inv} is the nominal (pseudo)inverse of {dir}, with {n+1}
  rows and {m+1} columns. Except for roundoff errors, and a uniform
  positive scaling, the matrices satisfy {dir*inv*dir = dir} and
  {inv*dir*inv = inv}. */

hrn_pmap_t hrn_pmap_alloc(uint32_t m, uint32_t n);
  /* Allocates a projetive function from {T^m} to {T^n},
    initializes it with the (pseudo)identity. */

void hrn_pmap_free(hrn_pmap_t M);
  /* De-allocates the arrays of projective map {M}. */

void hrn_map_point(double p[], hrn_pmap_t *M, double q[]);
  /* Applies the projective function {M} to point {p[0..m]} of {T^m} 
   and stores the result (a point of {T^n}) into {q[0..n]};
   where {m = M.m} and {n = M.n}. */

void hrn_map_hyperplane(double h[], hrn_pmap_t *M);
  /* Applies the projective function {M} to the hyperplane {h[0..m]} of {T^m}. */

void hrn_pmap_compose(hrn_pmap_t *M, hrn_pmap_t *N, hrn_pmap_t *R);
  /* Sets {R} to the composition of {M} and {N}, applied in that order.
    The client must allocate {R} before calling this procedure. */

hrn_pmap_t hrn_pmap_inv(hrn_pmap_t *M);
  /* Returns the inverse of the projective function {M}.
    Namely swaps {M.dir} with {M.inv}, and {M.m} with {M.n}.
    WARNING: the result shares with {M} the arrays {dir} and {ind}. */
    
/* CANONICAL SIMPLEX */

void hrn_canonical_simplex(uint32_t d, uint32_t n, double p[]);
  /* Stores into {p[0..(d+1)*(n+1)-1}]} the canonical {d}-dimensional
    simplex of {T^n}, for any {d <= n}. 
    
    Corner 0 of the simplex is the hither origin {[1,0,0,..0]}, and
    corners 1 to {d} are the points at infinity of the first {d}
    coordinate axes. Thus, if {d < n} the simplex lies in the
    canonical embedding of {T^d} in {T^n}.
    
    The array {p} is interpreted as a {d+1} by {n+1} matrix, where
    each row is the homogeneous cooordinate vector of a corner of the
    simplex. So, this procedure actually stores into {p} the first
    {d+1} rows of the identity matrix.
    
    Alternatively, row {i} of the matrix can be interpreted as the
    homogeneous coefficient vector of the hyperplane that supports the
    facet of the simplex that is opposite to corner {i}. If {d < n},
    that hyperplane is parallel to the last {n-d} coordinate axes. */

/* REGULAR SIMPLEX */

void hrn_regular_simplex(uint32_t n, double p[]);
  /* Stores into {p[0..(n+1)^2-1}]} the geometry of a regular
    {n}-dimensional simplex of {T^n}. 
    
    The geometry of the simplex is that defined by
    {rmxn_regular_simplex} in {rmxn.h}, assuming the standard
    embedding of {R^n} as the hither part of {T^n}. Namely, the
    corners of this simplex are all finite and in the hither part, the
    circumradius is {sqrt(n)}, and the barycenter is at the origin
    {[1,0,0,...,0]}.   Corner number 0 is point {(-1,-1,..-1)} of
    {R^n}, while corner number {i}, for {i} in {1..n}, 
    {(1+n*c)u_{i-1}-(c,c,..c)}; where {u_i} is the unit point of Cartesian
    axis {i}, and {c = (sqrt(n+1)-1)/n}.
    
    The procedure expects {p} to be the address of an {n+1} by {n+1}
    matrix, and puts into row {i} the homogeneous coordinate vector of
    corner {i} of the simplex. The first coordinate (homogeneous
    weight) of all corners is 1.
    
    Alternatively, row {i} of the matrix {p} can be viewed as the
    homogeneous coefficients of the hyperplane that supports the facet
    opposite to corner {i}. More precisley, the dot product (in
    {R^{n+1}}) of row {i} and row {j} is 0 if {i != j}, and {n+1} if
    {i = j}. In other words, the product of the matrix by its
    transpose is {n+1} times the identity matrix.
    
    Let {x} be a finite point of {T^n} with homogeneous coordinates
    {x[0..n]}. The {n+1} barycentric coordinates {b[0..n]} of {x}
    relative to this simplex can be computed as {(p*x)/(n+1)/|x[0]|}.
    Note that these barycentric coordinates add to {+1} for hither
    points, and to {-1} for yonder points. */

#endif
