#ifndef aa_trapez_H
#define aa_trapez_H

/* Conversion of affine arithmetic forms to trapezoids. */
/* Last edited on 2008-01-18 21:43:46 by stolfi */

#include <ia.h>
#include <ia_trapez.h>
#include <aa.h>

ia_trapez_t aa_trapez_from_pair(Interval *xr, AAP xf, AAP yf);
  /* Computes an {ia_tarpez_t} that encloses the joint range of two AA
    forms {xf,yf}, clipped to the vertical band whose X projection is
    the interval {xr}. The centroid of {xf} must be inside {xr}. */

#endif
