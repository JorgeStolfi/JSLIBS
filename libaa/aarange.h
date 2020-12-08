/* Computes joint ranges of two or more affine forms. */
/* Last edited on 2003-07-28 02:47:40 by stolfi */

#ifndef aarange_H
#define aarange_H

#include <aa.h>
#include <flt.h>

void aa_2d_range
  ( AAP x, 
    AAP y, 
    AATermCount *nv,
    double *xv,
    double *yv
  );
  /* Computes the vertices of the joint range of {x} and {y}.
    Returns the number of vertices in {nv}, and the vertices in
    {xv[0..2*nv-1],yv[0..2*nv-1]}.  These vectors should be 
    allocated with at least {(x->np)+(y->np)} elements.
    If the joint range degenerates to a single point, 
    returns {*nv = 0} and stores nothing in {xv,yv}. */

#endif
