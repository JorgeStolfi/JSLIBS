#ifndef psp_basic_H
#define psp_basic_H

/* Basic definitions for polynomial spline pulses on generic meshes. */
/* Last edited on 2011-09-19 00:57:43 by stolfilocal */

#define _GNU_SOURCE
#include <stdint.h>

#include <spline_basic.h>
#include <rgrid_basic.h>
#include <bz_basic.h>

/*
  REGULAR UNIDIMENSIONAL POLYNOMIAL SPLINES 
  
  A /unidimensional regular polynomial spline/ (RUPS) is a {P}-spline over a 
  regular unidimensional grid {G}. */

typedef spline_degree_t psp_degree_t;
  /* Degree of a unidimensional polynomial spline. */

/* 
  CONTINUITY OF A UNIDIMENSIONAL POLYNOMIAL SPLINE 
  
  A regular unidimensionals polynomial spline is continuous to all 
  orders inside each cell.  At each vertex, it may be continuous
  to all orders, or only to some finite order.  
 
  CONTINUOUS SPLINE SUBSPACES
    
  If {c} is any integer not less than {-1}, we denote by {S_c(F,G)} the space of
  all splines in {S(F,G)} that are continuous at least to order {c} with respect
  to all directional derivative operators. In particular, {S_0(F,G)} is the
  space of all continuous and piecewise smooth splines in {S(F,G)}.
 
  The {k}th derivative of a spline with continuity order {c} is defined
  everywhere if {k <= c}; for {k > c}, it is defined only inside the cells, not
  at the grid vertices. In particular, a spline with continuity {c = -1} is
  undefined at the grid vertices.
  
  If the domain {D} is a subset of {R^d} and {c} is an integer vector
  with {d} elements, then {S_c(F,G)} denotes the space of all
  {F}-splines on {G} that are continuous of order {c[i]} under
  differentiation with respect to domain coordinate number {i}.

*/
  
typedef spline_cont_t psp_cont_t;
  /* Continuity order of a polynomial spline. */
  
/*
  MULTI-DIMENSIONAL POLYNOMIALS
  
  Let {g[0..d-1]} be a vector of {d} natural numbers. A /polynomial of
  degree {g}/ is a function {f} from {R^d} to {R} such that, for any {i} in
  {0..d-1} and any {x} in {R^d}, the function {t --> f(x + t*e^d_i)} 
  is a polynomial of degree {g[i]} on the parameter {t}.
  
  We denote by {P^g} the vector space of all such polynomials.
  
  POLYNOMIAL SPLINES
  
  A /polynomial spline of degree {g}/ is a spline on the family {P^g},
  i.e. an element of {S(P^g,G)} for some mesh {G}. */

/*
  
  MULTI-DIMENSIONAL SPLINES
  
  In what follows, we consider the special case of meshes whose domain
  {D} is a subset of {R^d}, for some positive integer {d}, with the
  natural topology of the latter.
  
  For any {d} and any {i} in {0..d-1}, let {e^d_i} denote the unit
  vector of coordinate axis number {i} in {R^d}.
  
  POLYNOMIAL SPLINES ON REGULAR GRIDS
  
  A /polynomial spline/ on a regular {d}-dimensional mesh {G} with domain {\RR^d} 
  is a polynomial function from {\RR^d} 
  
  If {f} is a polynomial and {x} is a point of a cycle {D} with period {h}, we
  will interpret {f(x)} as equal to {f(\eta_h(x))}.
  
  A /unidimensional regular polynomial spline/ (RUPS) is a {P}-spline over a 
  regular unidimensional grid {G}.
  
  CONTINUOUS POLYNOMIAL SPLINE SPACES

  If {f} is defined on a {d}-dimensional domain, and {c[0..d-1]} is a vector
  of {d} integers, then saying that {f} is continuous of multi-order {c} means
  that, for almost any {x}, the univariate function {g(t) = f(x+t*u[i])} is
  continuous of order {c[i]} for any {i} in {0..d-1}, where {u[i]} is the
  direction vector of cordinate axis {i}.

  We are particularly interested in the space {S_c(P^g,G)} 
  for some given mesh {G}; namely, the polynomial splines on {G}
  with pieces of degree {g[0..d-1]} with continuity {c[0..d-1]}. */

#endif
