#ifndef minn_quad_H
#define minn_quad_H

/* Tools for optimizing {d}-dimensional function by enumeration. */
/* Last edited on 2023-03-27 15:09:28 by stolfi */ 

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <minn.h>

void minn_quad
  ( int32_t n,        /* Dimension of search space. */
    minn_goal_t *F,   /* Function to be minimized. */
    double dMax,      /* Radius of search domain, or {+INF}. */
    bool_t box,       /* True to search in the cube, false in the ball. */
    double tol,       /* Desired precision. */
    double v[],       /* (OUT) Minimum vector found. */
    double *Fval_P    /* (OUT) Goal function value at the minimum. */
  );
  /* Finds a vector {v[0..n-1]} of {\RR^n} that minimizes the goal function {F(n,v)}
    over a search domain {\RD} defined by {dMax} and {box}.
    
    Uses a quadratic approximation method of optimization, specifically
    the iterated SVE (Simplex Vertex-Edge) with some max number of iterations. 
    
    If {dMax} is {+INF}, the search domain {\RD} is the whole of {\RR^n},
    and {box} is ignored. Otherwise, if {box} is true the search domain
    is the signed cube {[-dMax _ +dMax]^n}, else it is the ball 
    {{v \in \RR^n : |v| <= dMax}}.
    
    The parameter {tol} must be a positive number. The procedure will attempt 
    find an approximation to the minimum with Euclidean error {tol} or less.
    
    On output, {v[0..n-1]} will be the best vector found. The value
    of {F(n,v)} is returned on {*Fval_P}. */

#endif
