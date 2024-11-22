/* qmin_simplex.h - quadratic minimization in a simplex. */
/* Last edited on 2024-11-20 08:58:50 by stolfi */

#ifndef qmin_simplex_H
#define qmin_simplex_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

/* In all the procedures below, two-dimensional matrices are stored
  into one-dimensional vectors, in row-by-row order. That is, an {m×n}
  matrix {A[0..m-1,0..n-1]} is stored as a vector {A[0..m*n-1]}, 
  with entry {A[i,j]} of the matrix in element {A[n*i+j]} of the
  vector. */
    

void qms_quadratic_min(uint32_t n, double A[], double b[], double x[]);
  /* Finds the minimum argument {x} of a quadratic function 
    {Q(x) = x' A x - 2 x'b + c}, subject to the constraints
    {x[i] >= 0} for all {i}; where {A} is a known positive
    semidefinite {n × n} matrix, {b} is a known {n}-vector, {x} is an
    unknown {n}-vector, and {x'} denotes the transpose of {x}.

    The optimum argument is returned in the vector {x}, whose original
    contents is ignored. The arrays {A} and {b} are not modified.
    
    Observe that the gradient of {Q} with respect to {x} is {A x - b}.
    It follows that the computed solution vector {x} is defined by the
    following conditions, for all {i} in {0..n-1}:
      
      (1): {x[i] >= 0} for all {i}.
      (2): if {x[i] > 0}, then {(A x - b)[i] ~ 0.0}.
      (3): if {x[i] == 0}, then {(A x - b)[i] >= 0.0}.
    
    Condition (2) means that the gradient must be zero along axis
    {i} when {x[i]} is far from the constraint {x[i] >= 0}; and
    condition (3) says that the gradient must be directed away from
    that constraint when {x[i]} is blocked by it. */

#endif
