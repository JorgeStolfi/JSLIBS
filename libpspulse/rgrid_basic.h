#ifndef rgrid_basic_H
#define rgrid_basic_H

/* Regular grids, infinite or toroidal. */
/* Last edited on 2020-10-03 21:48:22 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <spline_basic.h>

#include <interval.h>
#include <box.h>
#include <bool.h>
#include <sign.h>
#include <vec.h>
#include <indexing.h>

/* 
  SIMPLE UNIDIMENSIONAL DOMAINS
  
  A /cycle/ or /flat circle/ is an interval {[0_omega)}, for some positive real
  number {omega} (the /period/ or /extent/ of the cycle), endowed with the topology
  that results from assuming that the superior endpoint {omega} is identical to the
  inferior endpoint {0}.  By definition it is also a metric space, with the 
  metric 
  
    {dist(x,y) = MIN{ |x - y|, |x - y - omega|, |x - y + omega| }}
  
  The flat circle {[0_omega)} is homomorphic to the set of reals modulo the
   equivalence relation that equates any two reals which differ by an
  integer multiple of {omega}.
  
  A /simple unidimensional domain/ (SUD) is either the real line or a cycle. The
  /period/ of the real line, by definition, is {+oo}.
  
  The /canonical embedding/ {\eta_h} of any cycle {D} with period {h}   
  maps each point {x} of {D} to the unique real number {z} such that
  {0 <= z < h} and {x = { z + k*h : k in \RZ }}.  Note that {\eta_h} is
  one-to-one, and isometric except when {x} is a multiple of {h}.
  
  REGULAR UNIDIMENSIONAL GRIDS
  
  A /regular unidimensional grid/ (RUG) is a partition of a simple
  unidimensional domain {D} into a non-empty set of singletons, including {{0}},
  and a set of finite open intervals of equal length. 
  
  These subsets are the /items/ of the grid; the singletons are its /vertices/ or
  /knots/ and the intervals are its /cells/.
  
  In other words, a RUG is a mesh on a SUD (as defined in {spline_basic.h})
  that has {{0}} as a vertex and congruent cells. 
  
  The common finite length {delta} of the cells is the /step/ of the grid.
  
  Note that, in an RUG, all items with the same dimension are congruent and have
  topologically similar contexts. In particular, every vertex (including the
  origin 0) separates two cells, and every cell (including the last cell
  {[wd-delta _ wd)}) is bounded by two vertices. (However, if the grid has a
  single cell, then the two ends of that cell are the same vertex, and that
  vertex has the same cell on both sides.)
  
  SIZE OF A UNIDIMENSIONAL GRID
  
  The /size/ of the grid is the number {gsz = omega/delta} of cells. Note 
  that {gsz} is a positive integer if the domain is a cycle, or {+oo} 
  if it is the real line.  */
   
typedef ix_size_t rgrid_size_t;
  /* Size of a regular unidimensional grid, in multiples of the grid step.
    By convention, the value 0 denotes an infinite size. */
 
#define rgrid_size_MAX (ix_MAX_SIZE)
  /* Max size of a grid, to simplify avoidance of overflow. */

#define rgrid_size_FMT "%lu"
  /* Suitable format to print an {rgrid_size_t} value. */

/*
  POSITIONS IN A UNIDIMENSIONAL GRID
  
  A vertex {v} on an unidimensional grid can be identified by an
  integer /position/ {pos}, the displacement from the origin to {p}
  divided by the corresponding grid step {delta[i]}.  In particular, the
  vertex {{0}} has position 0.
  
  Similarly, a cell can be identified by the position of its inferior
  endpoint. In particular, the cell whose inferior endpoint is 0 (the
  /base cell/ of the grid) has position 0.
    
  Note that if the domain is periodic, positions that differ by
  a multiple of the grid size {gsz} denote the same item. */

typedef ix_pos_t rgrid_pos_t;
  /* Displacement of an item in a regular unidimensional grid, or a relative
    displacement between two items; both in multiples of the grid step. Also. */
 
#define rgrid_pos_MAX (ix_MAX_POS)
  /* Max absolute position in a grid, to simplify avoidance of overflow. */

#define rgrid_pos_FMT "%ld"
  /* Suitable format to print an {rgrid_pos_t} value. */

/*
  SIMPLE MULTIDIMENSIONAL DOMAINS
  
  A /simple multidimensional domain/ (SMD) is a topological space {D} that is
  the Cartesian product of a finite number {d} of simple unidimensional domains 
  {D[0],.. D[d-1]}, the /coordinate axes/ of {D}.  
  
  Such a space {D} is obviousy a {d}-dimensional manifold.  If all coordinate
  axes are cycles, then {D} is a {d}-dimensional torus.  If all are real lines,
  then {D} is simply {\RR^d}.  More generally, if {D[i]} is a cycle, we say that
  {D} has /toroidal/, /cyclic/, or /circular topology/ along axis {i}.
  
  The /period vector/ of {D} is the real vector {omega[0..d-1]} where {omega[i]} is the
  period of {D[i]}.

  REGULAR MULTIDIMENSIONAL GRIDS
  
  A /regular multidimensional grid/ (RMG) is essentially the Cartesian product
  of regular unidimensional grids.
  
  Namely, let {G[0],..G[d-1]} be a list of {d} RUGs. The RMG defined by them is
  a mesh {G}, such that evey item {X} of {G} is the cartesian product
  {X[0]×..X[d-1]} of {d} items, each {X[i]} being an item of {G[i]}. Note that
  the dimension {k} of this item is the number of factors {X[i]} that are cells.
  
  Note also that items of {G}, thus defined, are a partition of the SMD
  {D=D[0]×..D[d-1]}, where each {D[i]} is the domain of {G[i]}. Note also that
  {D} is a {d}-dimensional SMD.
  
  In particular, the vertices of {G} are cartesian products of {d} vertices,
  each from one {G[i]}. One of these vertices is the singleton containing the
  origin {(0,..0)} of {\RR^d}.  Likewise, the cells of {G} are {d}-dimensional
  boxes, each the Cartesian product of {d} cells, one from each {G[i]}.
  
  STEP AND SIZE VECTORS
  
  The /step/ of {G} is the real vector {delta[0..d-1]} where {delta[i]} is the
  step of {G[i]}.
  
  The /size/ of {G} is the vector {gsz[0..d-1]} where {gsz[i]} is the size of
  {G[i]} (either an integer or {+oo}, encoded as 0).
  
  GRID POSITION VECTORS
  
  A vertex {v} on a fixed grid can be identified by an integer
  /position vector/ {pos[0..d-1]}, where {pos[i]} is the displacement
  from the origin to {v} along the coordinate axis {i}, divided by the
  corresponding grid step {delta[i]}. In particular, the vertex at the
  origin of {R^d} has position vector {(0,..0)}.
  
  Similarly, a cell can be identified by the position of its inferior
  corner. In particular, the cell whose inferior corner is the origin
  (the /base cell/ of the grid) has position vector {(0,..0)}.
  
  SPANNING AND NORMAL AXES OF AN ITEM 
  
  Each {k}-dimensional item {X} of {G} is a {k}-dimensional open box in {R^d}.
  (see {box.h}). Therefore, the axis numbers {0..d-1} can be partitioned into a
  subset {Spn(X)} of {k} /spanning axes/ (parallel to the box) and a subset
  {Nrm(X)} of {d-k} /normal axes/ (orthogonal to the box).
  
  Either of these two complementary sets define the /orientation/ of the item in
  {d}-space. All items of {G} with the same orientation are congruent translates
  of each other. Indeed, any item of {G} can be identified by its set of
  spanning axes, and by the position vector of its inferior corner.
  
  Note also that every item with dimension {k < d} lies on the boundary of
  {2*(d-k)} items with dimension {k+1}. Dually,every item with dimension {k > 0}
  is bounded by {2*k} items with dimension {k-1}. */

typedef box_dim_t rgrid_dim_t;
  /* Dimension of a domain or range space, mesh element, etc.. */
  
typedef box_axis_t rgrid_axis_t; 
  /* Identifies a coordinate axis of a multidimensional 
    space, from 0 to {d-1} where {d} is the space's dimension. */

typedef box_axis_set_t rgrid_axis_set_t; 
  /* A packed subset of coordinate axes (each a {rgrid_axis_t}). */

typedef box_axis_index_t rgrid_axis_index_t; 
  /* Identifies a coordinate axis among a set of axes: 0 is the lowest. */

#define rgrid_axis_NONE (box_axis_NONE)
  /* A {rgrid_axis_t} value that means "no axis". */

#define rgrid_axis_NOT_FOUND (box_axis_NOT_FOUND)
  /* A {rgrid_axis_t} value that means "no such axis" */
  
rgrid_axis_set_t rgrid_axis_set_complement(rgrid_dim_t d, rgrid_axis_set_t A);
  /* The complement of {A} relative to the set {0..d-1}. */

#endif
