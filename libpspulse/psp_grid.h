#ifndef psp_grid_H
#define psp_grid_H

/* Regular grids, infinite or toroidal. */
/* Last edited on 2014-05-15 22:49:49 by stolfilocal */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <interval.h>
#include <box.h>
#include <rgrid_basic.h>
#include <bz_basic.h>
#include <bool.h>
#include <sign.h>
#include <vec.h>

typedef rgrid_dim_t psp_dim_t;
  /* Number of coordinate axes in a domain or range space, mesh element, etc.. */

typedef rgrid_size_t psp_grid_size_t;
  /* Size of a grid, in multiples of the grid step.
    By convention, the value 0 denotes an infinite grid. */
 
#define psp_grid_size_MAX rgrid_size_MAX
  /* Max size of a grid, to simplify avoidance of overflow. */

#define psp_grid_size_FMT rgrid_size_FMT
  /* Suitable format to print an {psp_grid_size_t} value. */

typedef rgrid_pos_t psp_grid_pos_t;
  /* Displacement of an item in a regular unidimensional grid, or a relative
    displacement between two items; both in multiples of the grid step. */
  
#define psp_grid_pos_MAX rgrid_pos_MAX
  /* Max absolute position in a grid, to simplify avoidance of overflow. */

#define psp_grid_pos_FMT rgrid_pos_FMT
  /* Suitable format to print an {psp_grid_pos_t} value. */

typedef rgrid_axis_t psp_axis_t; 
  /* Identifies a coordinate axis of a multidimensional 
    space, from 0 to {d-1} where {d} is the space's dimension. */

typedef rgrid_axis_set_t psp_axis_set_t; 
  /* A packed subset of coordinate axes (each a {psp_axis_t}). */

typedef rgrid_axis_index_t psp_axis_index_t; 
  /* Identifies a coordinate axis among a set of axes: 0 is the lowest. */

#define psp_axis_NONE (rgrid_axis_NONE)
  /* A {psp_axis_t} value that means "no axis". */

#define psp_axis_NOT_FOUND (rgrid_axis_NOT_FOUND)
  /* A {psp_axis_t} value that means "no such axis" */
  
psp_axis_set_t psp_axis_set_complement(psp_dim_t d, psp_axis_set_t A);
  /* The complement of {A} relative to the set {0..d-1}. */


#endif
