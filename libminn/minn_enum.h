#ifndef minn_enum_H
#define minn_enum_H

/* Tools for optimizing {d}-dimensional function by enumeration. */
/* Last edited on 2025-04-01 09:10:13 by stolfi */ 

#include <stdint.h>

#include <bool.h>
#include <minn.h>

void minn_enum
  ( uint32_t n,       /* Dimension of search space. */
    minn_goal_t *F,   /* Function to be minimized. */
    bool_t dBox,      /* True to search in the cube, false in the ball. */
    double tol[],     /* Desired precision along each axis. */
    double v[],       /* (OUT) Minimum vector found. */
    double *Fval_P    /* (OUT) Goal function value at the minimum. */
  );
  /* Finds a vector {v[0..n-1]} in a certain domain {\RD}
    that minimizes the goal function {F(n,v)}.
    Uses exhaustive enumeration of a regular grid vectors {u} in
    the domain {\RD} with spacing {tol[i]} along each axis {i}.
    
    If {dBox} is true, the domain {\RD} is the axis-aligned {n}-dimensional cube
    {[-1 _ +1]^n}. If {dBox} is false, the domain is the {n}-dimensional ball
    of unit radius.
    
    On output, {v[0..n-1]} will be the best vector found. The value
    of {F(n,v)} is returned on {*Fval_P}.
    
    The probe points {q} will comprise a compact regular lattice with
    closest distance {tol} that includes the origin. Currently the
    lattice is an orthogonal grid with axes aligned to the coordinate
    axes. */

#endif
