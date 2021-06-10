/* rmxn_extra.h --- additional operation on MxN matrices */
/* Last edited on 2021-06-09 20:33:25 by jstolfi */

#ifndef rmxn_extra_H
#define rmxn_extra_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <rn.h>
#include <rmxn.h>

void rmxn_perturb_unif(int32_t m, int32_t n, double mag, double M[]);
  /* Adds an independent random perturbation of with uniform distribution 
    in {[-mag _ +mag]} to each element of {M}, which is assumed to
    have {m} rows and {n} columns. */ 

void rmxn_throw_ortho(int32_t n, double M[]);
  /* Stores in {M} a random positive orthonormal {n×n} matrix. 
    Its rows and columns are pairwise orthogonal and with length 1, and ditto
    for its columns; and its determinant is {+1}.  Assumes {M} has {N^2} elements. */

void rmxn_spin_rows(int32_t m, int32_t n, double A[], double M[]);
  /* Applies a random rotation to each row of {A}, which is assumed to
    have {m} rows and {n} columns. Equivalent to computing {M=A*N}
    where {N} is a random orthonormal {n×n} matrix. */

void rmxn_spin_cols(int32_t m, int32_t n, double A[], double M[]);
  /* Applies a random rotation to each column {A}, which is assumed to
    have {m} rows and {n} columns. Equivalent to computing {M=N*A}
    where {N} is a random orthonormal {m×m} matrix. */

void rmxn_shift_rows(int32_t m, int32_t n, double A[], double v[], double M[]);
  /* Adds the vector {v[0..n-1]} to each row of the matrix {A},
    yielding matrix {M}. Both matrices are assumed to have {m} rows
    and {n} columns. */

void rmxn_shift_cols(int32_t m, int32_t n, double A[], double v[], double M[]);
  /* Adds the vector {v[0..m-1]} to each column of the matrix {A},
    yielding matrix {M}. Both matrices are assumed to have {m} rows
    and {n} columns. */

/* THE CANONICAL {d}-SIMPLEX */

void rmxn_canonical_simplex(int32_t d, int32_t n, double V[]);
  /* Stores into {V[0..(d+1)*n-1]} the coordinates of the vertices
    of the /canonical {d}-dimensional simplex of {R^n}/, for any {d < n}.

    The canonical {d}-simplex has {d+1} corners; corner number {i} 
    is the unit point {u_i} of axis {i} of {R^n}. The array
    {V} is interpreted as a {d+1} by {n} matrix, where each row is
    a corner of the simplex. Thus, this function merely stores into
    {V} the first {d+1} rows of the {n} by {n} identity matrix.
    
    The canonical {d}-simplex is a regular simplex, with edge length
    {sqrt(2)} and circumradius {sqrt((d+1)/d)}. Its center is the
    point {(c,c,...c)} where {c = 1/(d+1)}. It is actually contained
    in the {d}-dimensional affine subspace of {R^n} consisting of all
    points whose first {d+1} coordinates add to 1, and whose {n-d-1}
    coordinates are all zero. Its isometric symmetries are all the
    {(d+1)!} permutations of the first {d+1} axes of {R^n}. */

double rmxn_canonical_simplex_radius(int32_t d);
  /* Circum-radius {rd(d) = sqrt(d/(d+1)} of the canonical
    {d}-simplex. Namely, the distance from any unit point {v} of
    {R^{d+1}} to the barycenter of all its unit points INCLUDING {v}.
    In particular,
      
      | {d}     |         0 |         1 |         2 |         3 |
      +---------+-----------+-----------+-----------+-----------+
      | {or(d)} |         0 | sqrt(1/2) | sqrt(2/3) | sqrt(3/4) |  */

double rmxn_canonical_simplex_subradius(int32_t d, int32_t k);
  /* Radius {sr(d,k) = sqrt((d-k)/((d+1)*(k+1))} of the
    {(d-1)}-dimensional sphere that is concentric with the canonical
    {d}-simplex and contains the center of all its {k}-dimensional
    faces.
    
    In particular, {sr(d,0)} is the circum-radius {rd(d)}; {sr(d,1)}
    is the distance from the simplex center to the midpoint of each
    edge; {sr(d,d-1)} is the in-radius; and {sr(d,d)} is zero. */

double rmxn_canonical_simplex_edge(int32_t d);
  /* Edge length {ed(d)} of the canonical {d}-simplex, namely
    {sqrt(2)} (even if {d == 0}). */
  
double rmxn_canonical_simplex_height(int32_t d);
  /* Height {ht(d) = sqrt((d+1)/d)} of the canonical {d}-simplex.
    Namely, the distance from any unit point {v} of {R^{d+1}} to the
    barycenter of the remaining {d} unit points DISTINCT from {v}. In
    particular,
      
      | {d}     |         0 |         1 |         2 |         3 |
      +---------+-----------+-----------+-----------+-----------+
      | {ht(d)} |       +oo | sqrt(2/1) | sqrt(3/2) | sqrt(4/3) | */

double rmxn_canonical_simplex_measure(int32_t d);
  /* The {d}-dimensional measure {ms(d) = sqrt(d+1)/(d!)} of the
    canonical {d}-simplex. In particular, 
      
      | {d}    |         0 |         1 |         2 |         3 |
      +--------+-----------+-----------+-----------+-----------+
      | {ms(d) |         0 | sqrt(2)/1 | sqrt(3)/2 | sqrt(4)/6 | */
  
/* A REGULAR {n}-SIMPLEX */

void rmxn_regular_simplex(int32_t n, double V[]);
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

double rmxn_regular_simplex_radius(int32_t n);
  /* Circum-radius {RD(n) = sqrt(n)} of the regular {n}-simplex above. */

double rmxn_regular_simplex_subradius(int32_t n, int32_t k);
  /* Radius {SR(n,k) = sqrt((n-k)/(k+1))} of the {(n-1)}-dimensional
    sphere in {R^n} that is concentric with the canonical {n}-simplex
    above and contains the center of all its {k}-dimensional faces.
    
    In particular, {SR(n,0)} is the curcum-radius {RD(n)}; {SR(n,1)}
    is the distance from the simplex center to the midpoint of each
    edge; {SR(n,n-1)} is the in-radius; and {SR(n,n)} is zero. */

double rmxn_regular_simplex_edge(int32_t n);
  /* Edge length {ED(n) = sqrt(2*(d+1))} of the regular {n}-simplex above.
    Arbitrarily set to {sqrt(2)} if {d == 0}. */

double rmxn_regular_simplex_height(int32_t n);
  /* Height {HT(n) = (d+1)/sqrt(d)} of the regular {n}-simplex above. */
  
double rmxn_regular_simplex_measure(int32_t n);
  /* The {n}-dimensional measure {MS(n) = (d+1)^{(d+1)/2}/(d!)} of the
    regular {n}-simplex above. */

/* BARYCENTRIC COORDINATES */

void rmxn_throw_canonical_simplex(int32_t d, double x[]);
  /* Generates a random point {x[0..d]} uniformly distributed
    in the canonical {d}-dimensional simplex of {R^{d+1}}.
    Assumes that {x} has {d+1} elements. */
     
void rmxn_throw_canonical_simplex_ball(int32_t d, double x[]);
  /* Generates a random point {x[0..d]} in the {d}-dimensional 
    ball circumscribed on the {d}-dimensional canonical simplex.
    Assumes that {x} has {d+1} elements. */

#endif

