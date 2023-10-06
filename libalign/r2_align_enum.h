#ifndef r2_align_enum_H
#define r2_align_enum_H

/* Tools for optimizing translational alignment of 2D objects by enumeration. */
/* Last edited on 2023-09-07 18:05:53 by stolfi */ 

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <r2.h>
#include <i2.h>
#include <r2_align.h>

void r2_align_enum
  ( int32_t ni,              /* Number of objects to align. */
    r2_align_mismatch_t *F2, /* Function that evaluates the mismatch between the objects. */
    r2_t arad[],             /* Max delta vector coordinates for each object. */
    bool_t bal,              /* True if alignment vector adjustments should be balanced. */
    double tol,              /* Desired precision. */
    r2_t p[],                /* (IN/OUT) Initial and optimized alignment vector. */
    double *F2val_P          /* (OUT) Mismatch for the computed alignment vector. */
  );
  /* Adjusts an alignment vector {p[0..ni-1]} for some {ni} objects so
    as to minimize the mismatch function {F2(ni,p)}.
    
    On input, {p[0..ni-1]}, must be a guess for the optimum alignment.
    On output, {p[0..ni-1]} will be the best alignment found. The value
    of {F2(ni,p)} is returned on {*F2val_P}.
    
    Uses exhaustive enumeration of alignment vectors {q} in the search
    ellipsoid {\RF} derived from the basic ellipsoid {\RE} whose center
    is the initial value of {p} and whose radius vector is {arad}, with
    the balancing constraint if {bal} is true.
    
    The probe points {q} will comprise a compact regular lattice with
    closest distance {tol} that includes the initial guess {p}.
    Currently the lattice is an orthogonal grid with axes aligned to the
    major axes of the ellipsoid {\RF}. */

#endif
