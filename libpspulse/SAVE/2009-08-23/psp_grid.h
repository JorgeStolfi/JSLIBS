/* Regular grids, infinite or toroidal. */
/* Last edited on 2009-08-23 17:50:06 by stolfi */

#ifndef psp_grid_H
#define psp_grid_H

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <interval.h>
#include <bool.h>
#include <sign.h>
#include <vec.h>

/* 
  UNIFORM UNIDIMENSIONAL GRIDS
  
  A /uniform unidimensional grid/ is a partition {G} of the reals {R} into
  open intervals and single real numbers, the /items/ of the grid,
  such that (1) one of the items is the single number {0}, and (2) all
  intervals have the the same length {delta} -- the /step/ of the
  grid.
  
  The singleton items and the interval items are also called the
  /vertices/ and /cells/ of the grid, respectively.
  
  Note that every item has the same geometry and topology as any other
  item of the same dimension. In particular, every vertex is incident
  to two cells, and every cell is bounded by two vertices.
  
  GRID POSITIONS
  
  A vertex {v} on an unidimensional grid can be identified by an
  integer /position/ {pos}, the displacement from the origin to {p}
  divided by the corresponding grid step {delta[i]}.  In particular, the
  vertex at 0 has position 0.
  
  Similarly, a cell can be identified by the position of its inferior
  endpoint. In particular, the cell whose inferior endpoint is 0 (the
  /base cell/ of the grid) has position 0. */

typedef int64_t psp_grid_pos_t;
  /* Displacement of a grid item, in multiples of the corresponding
    grid step. Also a relative displacement between two items. */
  
/*
  CIRCULAR TOPOLOGY
  
  A /flat circle/ is an interval {[0_wd)}, for some positive real
  number {wd} (the /period/ or /extent/ of the flat circle), endowed
  with the topology that results from assuming that the superior endpoint
  {wd} is identical to the inferior endpoint {0}.
  
  It can be understood as the quotient of the real line {R} by an
  equivalence relation that equates any two reals which differ by an
  integer multiple of {wd}.
 
  CIRCULAR GRIDS
  
  A /circular grid/ is a grid on a flat circle whose period {wd} is
  some integer multiple {gsz*delta} of the grid's step {delta}. The
  integer {gsz} is called the grid's /size/.
  
  Thus a circular grid consists of a finite set {G} of items, taken
  from some (infinite) unidimensional grid, whose union {D = \U G} is
  the interval {[0 _ wd)} The interval {D} is called the /domain/ of
  the grid.
  
  Note that in a circular grid {G}, as in an infinite grid, all items
  of the same dimension are topologically alike. Namely, every vertex
  (including the origin 0) separates two cells, and every cell
  (including the last cell {[wd-delta _ wd)}) is bounded by two
  vertices. (However, if the grid has a single cell, then the two ends
  of the cell are the same vertex, and that vertex has the same cell
  on both sides.) */
  
typedef uint64_t psp_grid_size_t;
  /* Size of a circular grid, in multiples of the grid step. */

/* 
  UNIFORM MULTIDIMENSIONAL GRIDS
  
  A /uniform multidimensional grid/ is essentially the Cartesian product of 
  unidimensional grids. */

/*
  MULTI-DIMENSIONAL GRIDS
  
  A /{d}-dimensional grid/ can be seen as the Cartesian product of {d}
  partitions that are unidimensional grids (see {udg_grid.h}). It is a
  partition {G} of {R^d} into boxes of various dimensions, the /items/
  of the grid, such that 
  
    (1) the intersection of the closure of any two items of {G} 
    is the union of items of {G};
    
    (2) for any axis {i}, the projections of all items parallel
    to {i} are open intervals with the same length {delta[i]}
    -- the /step/ of the grid along that  axis; and 
    
    (3) one of the items consists of a single point, the
    origin of {R^d}.
  
  A /cell/ of a grid is an item {C} with maximum dimension, i.e.
  {dim(C) == d}. Grid items with dimensions 0, 1, 2, 3 are called
  /vertices/,  /edges/, /faces/, and /blocks/ of the grid, respectively.
  
  It follows from this definition that all items with the same
  oriantation (same set of normal axes) are congruent translates of each
  other.  Indeed, any item {I} can be identified by its orientation and 
  the unique vertex that is the inferior corner of {I} (or, equivalently, by the 
  unique cell which has {I} as a inferior face).
  
  It folows also that every item with dimension {m < d} lies on the
  boundary of {2^{d-m}} items with dimension {m+1}.
  
  GRID POSITION VECTORS
  
  A vertex {v} on a fixed grid can be identified by an integer
  /position vector/ {pos[0..d-1]}, where {pos[i]} is the displacement
  from the origin to {v} along the coordinate axis {i}, divided by the
  corresponding grid step {delta[i]}. In particular, the vertex at the
  origin of {R^d} has position vector {(0,..0)}.
  
  Similarly, a cell can be identified by the position of its inferior
  corner. In particular, the cell whose inferior corner is the origin
  (the /base cell/ of the grid) has position vector {(0,..0)}.
  
  More generally, an item can be identified by its set of normal or
  spanning vectors, and by the position vector of its inferior corner. */
  
/*
  TOROIDAL TOPOLOGY
  
  A (/canonical/) /{d}-torus/ is the Cartesian product of {d} flat
  circles (see {cp_pulse.h}). It is the Cartesian product of {d}
  intervals of the form {[0_wd[i])}, which is a {d}-dimensional box
  with inferior corner {(0,0,...)} plus all the inferior faces of that box.
  The topology is such that every point on any superior face of that box
  is identified with the corresponding point on the opposite inferior face.
  
  Each number {wd[i]} is called the /period/ or /extent/ of the torus
  along axis {i}.

  It can be understood as the quotient of {R^d} by an equivalence
  relation that equates any two points of {R^d} which, along every
  coordinate axis {i}, differ by an integer multiple of {wd[i]}.
 
  TOROIDAL GRIDS
  
  A /{d}-dimensional toroidal grid/ can be seen as the Cartesian
  product (in the partition sense) of {d} unidimensional circular
  grids (see {cp_pulse.h}). It consists of a finite subset {G} of
  items of an (infinite) {d}-dimensional grid, whose union is a
  {d}-dimensional torus {D} (the the /domain/ of the grid).
  
  The period {wd[i]} of the torus along axis {i} must be some integer
  multiple {gsz[i]*delta[i]} of the grid's step {delta[i]} along that
  axis. The integer {gsz[i]} is called the grid's /size/ along axis
  {i}. When {d==1}, we often say /circular/ instead of toroidal.
  
  Note that in a toroidal grid, as in an infinite grid, all items of a
  given dimension and orientation are topologically equivalent. */

#endif
