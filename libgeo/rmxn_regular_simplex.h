/* rmxn_regular_simplex.h --- additional operation on MxN matrices */
/* Last edited on 2024-12-05 10:28:33 by stolfi */

#ifndef rmxn_regular_simplex_H
#define rmxn_regular_simplex_H

#include <stdio.h>
#include <stdint.h>

#include <rn.h>
#include <rmxn.h>

/* REGULAR {n}-SIMPLICES IN {n}-SPACE */

void rmxn_regular_simplex(uint32_t n, double V[]);
  /* Stores into {V[0..(n+1)*n-1]} the coordinates of a regular
    {n}-dimensional simplex in {R^n}/.

    This simplex has {n+1} corners. Corner 0 is the point
    {(-1,-1,..-1)} of {R^n}.  Corner number {i}, for {i} in {1..n}, is
    {(1+n*c)u_{i-1}-(c,c,..c)}; where {u_i} is the unit point of Cartesian
    axis {i}, and {c = (sqrt(n+1)-1)/n}. 
    
    The simplex is regular and centered at the origin. It has
    {(n+1)!} linear symmetries (isometric linear maps of of {R^n} that
    take the simplex to itself). In particular, the {n!} symmetries
    that leave vertex 0 fixed are the permutations of the coordinate
    axes of {R^n}. The remaining {(n+1)! - n!} symmetries are more
    complicated.
    
    The array {V} is interpreted as a matrix with {n+1} rows and {n}
    columns, where each row is a corner of the simplex.
    
    Each row {i} of {V} can be interpreted also as a vector of {R^n}
    that is perpendicular to the facet of the simplex that is opposite
    to corner {i}. The length of the vector is such that a point {x}
    of {R^n} is coplanar with that facet iff the scalar product of row
    {i} of {V} and {x} is {-1}.
    
    Let {x} be a finite point of {R^n} with coordinates {x[0..n-1]}.
    The {n+1} barycentric coordinates {b[0..n]} of {x} relative to
    this simplex can be computed as {((1,1,..1) + V*x)/(n+1)}. Note
    that these barycentric coordinates add to {+1}. */

double rmxn_regular_simplex_radius(uint32_t n);
  /* Circum-radius {RD(n) = sqrt(n)} of the regular {n}-simplex above. */

double rmxn_regular_simplex_subradius(uint32_t n, uint32_t k);
  /* Radius {SR(n,k) = sqrt((n-k)/(k+1))} of the {(n-1)}-dimensional
    sphere in {R^n} that is concentric with the canonical {n}-simplex
    above and contains the center of all its {k}-dimensional faces.
    
    In particular, {SR(n,0)} is the curcum-radius {RD(n)}; {SR(n,1)}
    is the distance from the simplex center to the midpoint of each
    edge; {SR(n,n-1)} is the in-radius; and {SR(n,n)} is zero. */

double rmxn_regular_simplex_edge(uint32_t n);
  /* Edge length {ED(n) = sqrt(2*(d+1))} of the regular {n}-simplex above.
    Arbitrarily set to {sqrt(2)} if {d == 0}. */

double rmxn_regular_simplex_height(uint32_t n);
  /* Height {HT(n) = (d+1)/sqrt(d)} of the regular {n}-simplex above. */
  
double rmxn_regular_simplex_measure(uint32_t n);
  /* The {n}-dimensional measure {MS(n) = (d+1)^{(d+1)/2}/(d!)} of the
    regular {n}-simplex above. */

#endif

