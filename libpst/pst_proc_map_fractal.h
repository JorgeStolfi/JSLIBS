#ifndef pst_proc_map_fractal_H
#define pst_proc_map_fractal_H

/* pst_proc_map_fractal.h -- procedures for creating procedurally-defined fractal images. */
/* Last edited on 2024-12-22 12:06:38 by stolfi */

#include <stdint.h>

#include <r2.h>
#include <pst_proc_map.h>

void pst_proc_map_fractal_cone(r2_t *p, double eps, double seed, double *z, r2_t *dz);
  /* Height map of a fractal mountain whose base is approximately the
    unit disc of the XY plane, with height 1 at the origin and
    dropping irregularly towards 0 along the base. 
    
    The function is actually a digital terrain with a semiregular (4,8) 
    mesh of isosceles right triangles, whose edges are smaller than {eps}.
    
    The {seed} parameter should be a number in {(0 _ 1)} used as a seed
    for pseudo-random number generators. */

#endif
