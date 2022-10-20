#ifndef r2_align_enum_H
#define r2_align_enum_H

/* Tools for optimizing translational alignment of 2D objects by enumeration. */
/* Last edited on 2022-10-20 06:16:52 by stolfi */ 

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <r2.h>
#include <i2.h>
#include <r2_align.h>

void r2_align_enum
  ( int32_t ni,                  /* Number of objects to align. */
    r2_align_mismatch_t *F2, /* Function that evaluates the mismatch between the objects. */
    r2_t arad[],             /* Max alignment adjustment for each object. */
    double tol,              /* Desired precision. */
    r2_t p[],                /* (IN/OUT) Corresponding points in each object. */
    double *F2valP           /* (OUT) Mismatch for the computed alignment vector. */
  );
  /* Adjusts an alignment vector {p[0..ni-1]} for some {ni} objects so
    as to minimize the mismatch function {F2(ni,p)}.
    
    On input, {p[0..ni-1]}, must be a guess for the optimum alignment.
    On output, {p[0..ni-1]} will be the best alignment found. The value
    of {F2(ni,p)} is returned on {*F2valP}.
    
    Uses exhaustive enumeration of alignment vectors {q} in the search ellipsoid 
    {\RF} derived from the basic ellipsoid {\RE} hose center is the initial value of
    {p} and whose radius vector is {arad}. 
    
    The probe points {q} will comprise a compact regular lattice with
    closest distance {tol} that includes the initial guess {p}.
    Currently the lattice is an orthogonal grid with axes aligned to the
    major axes of the ellipsoid {\RF}. */

#endif
