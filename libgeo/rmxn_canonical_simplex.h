/* rmxn_extra.h --- additional operation on MxN matrices */
/* Last edited on 2024-11-20 13:04:04 by stolfi */

#ifndef rmxn_canonical_simplex_H
#define rmxn_canonical_simplex_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <rn.h>
#include <rmxn.h>

/* THE CANONICAL {d}-SIMPLEX IN {d+1}-SPACE */

void rmxn_canonical_simplex(uint32_t d, uint32_t n, double V[]);
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

double rmxn_canonical_simplex_radius(uint32_t d);
  /* Circum-radius {rd(d) = sqrt(d/(d+1)} of the canonical
    {d}-simplex. Namely, the distance from any unit point {v} of
    {R^{d+1}} to the barycenter of all its unit points INCLUDING {v}.
    In particular,
      
      | {d}     |         0 |         1 |         2 |         3 |
      +---------+-----------+-----------+-----------+-----------+
      | {or(d)} |         0 | sqrt(1/2) | sqrt(2/3) | sqrt(3/4) |  */

double rmxn_canonical_simplex_subradius(uint32_t d, uint32_t k);
  /* Radius {sr(d,k) = sqrt((d-k)/((d+1)*(k+1))} of the
    {(d-1)}-dimensional sphere that is concentric with the canonical
    {d}-simplex and contains the center of all its {k}-dimensional
    faces.
    
    In particular, {sr(d,0)} is the circum-radius {rd(d)}; {sr(d,1)}
    is the distance from the simplex center to the midpoint of each
    edge; {sr(d,d-1)} is the in-radius; and {sr(d,d)} is zero. */

double rmxn_canonical_simplex_edge(uint32_t d);
  /* Edge length {ed(d)} of the canonical {d}-simplex, namely
    {sqrt(2)} (even if {d == 0}). */
  
double rmxn_canonical_simplex_height(uint32_t d);
  /* Height {ht(d) = sqrt((d+1)/d)} of the canonical {d}-simplex.
    Namely, the distance from any unit point {v} of {R^{d+1}} to the
    barycenter of the remaining {d} unit points DISTINCT from {v}. In
    particular,
      
      | {d}     |         0 |         1 |         2 |         3 |
      +---------+-----------+-----------+-----------+-----------+
      | {ht(d)} |       +oo | sqrt(2/1) | sqrt(3/2) | sqrt(4/3) | */

double rmxn_canonical_simplex_measure(uint32_t d);
  /* The {d}-dimensional measure {ms(d) = sqrt(d+1)/(d!)} of the
    canonical {d}-simplex. In particular, 
      
      | {d}    |         0 |         1 |         2 |         3 |
      +--------+-----------+-----------+-----------+-----------+
      | {ms(d) |         0 | sqrt(2)/1 | sqrt(3)/2 | sqrt(4)/6 | */
  
void rmxn_canonical_simplex_throw(uint32_t d, double x[]);
  /* Generates a random point {x[0..d]} uniformly distributed
    in the canonical {d}-dimensional simplex of {R^{d+1}}.
    Assumes that {x} has {d+1} elements. */
     
void rmxn_canonical_simplex_ball_throw(uint32_t d, double x[]);
  /* Generates a random point {x[0..d]} in the {d}-dimensional 
    ball circumscribed on the {d}-dimensional canonical simplex.
    Assumes that {x} has {d+1} elements. */

#endif

