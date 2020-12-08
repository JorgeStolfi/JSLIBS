#ifndef pst_geom_H
#define pst_geom_H

/* pst_geom.h -- misc geometric tools. */
/* Last edited on 2009-02-28 21:03:53 by stolfi */

#include <r2.h>
#include <r3.h>

#include <pst_basic.h>

/* UTILITIES */

void pst_geom_clip_dir(r3_t *udir, r3_t *sdir, double ard);
  /* Clips the direction vector {udir} to a spherical cap with center
    on the direction vector {sdir} and angular radius {ard. Namely, if
    the angle between {udir} and {sdir} is greater than {ard} radians,
    sets {udir} to the unit vector that lies on the shortest arc
    between the two, at {ard} radians from {sdir}. */

#endif
