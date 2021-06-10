/* Last edited on 2021-06-09 20:20:53 by jstolfi */

#include <stdint.h>

/* NON-CANONICAL REGULAR {n}-SIMPLICES */
  
void rmxn_regular_simplex(int32_t n, double V[], double U[]);
  /* This simplex also has {n+1} corners, and can be viewed as sitting
    either (1) in the finite part of the {n}-dimensional projective space
    {T^n}, or (2) in the {n}-dimensional subspace of Cartesian space
    {R^{n+1}} consisting of the points {x} with {x[0] == 1}; or, by
    discarding the first coordinate of every corner, (3) in Cartesian
    space {R^n}.
    
    In either of these interpretations, the result is isometric to the
    canonical {n}-simplex; hence it is regular, with edge length 2 and
    circumradius {n/sqrt(n+1)}.  In interpretation (3), its first {n} corners
    are {u_i + (a,a,...,a)}, and the last corner is {(b,b,...b)};
    where {b} negative.
    
    
  */
  
void rmxn_barycentric_matrix(int32_t n, double V[], double U[]);
  /* Stores into {U[0..(n+1)*(n+1)-1]} a projetive matrix that can be
    used to compute the barycentric coordinates of any given point of
    {R^n} relative to a given {n}-simplex {V}.
 
    The given array {V} is interpreted as an {(n+1)×n} matrix, where
    each row is a simplex corner. The returned array {U} is
    interpreted as an {(n+1)×(n+1)} projective matrix, where the first
    row is the translation term, and the last {n} rows and columns are
    an orthonormal transformation of {R^n}.
    
    After the call, the barycentric coordinates of a point {p} of
    {R^n} with respect to the simplex {V} can be obtained by the
    vector-matrix product {x*U}, where {x} is the homogeneous
    coordinate vector of {p} (with {x[0]} being the weight
    coordinate).
    
    Thus, {U*V} will be the identity matrix of {R^n},
    and {V*U} will be the orthogonal projection of {R^{n+1}} on the
    affine subspace spanned by the canonical {n}-simplex (the set of
    all points whose coordinates add to 1). */

