#ifndef psp_basic_H
#define psp_basic_H

/* Basic definitions for polynomial splines on generic meshes. */
/* Last edited on 2009-08-23 20:32:27 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <box.h>

/*
  GENERAL MESHES
  
  A /mesh/ is a finite partition {G} of some topological space {D} (the
  /domain/) such that (1) each element of {G} is a topological open
  {k}-ball of {D}, for some {k}; and (2) the topological closure in
  {D} of any element of {G} is the union of finitely many elements of {G}. 
  
  The elements with maximum dimension are the /cells/ of the mesh.
  Elements with dimension {k = 0,1,2,3} are called /vertices/, /edges/, 
  /panes/, and /blocks/, respectively.
  
  A /submesh/ is any subset of a mesh; the /domain/ of the submesh is 
  is the union of its elements. */
  
typedef box_dim_t psp_dim_t;
  /* Dimension of a domain or range space, mesh element, etc.. */
  
/*  
  GENERAL SPLINES
  
  Let {G} be a submesh with domain {D}, and {F} be a family of
  functions from {D} to some vector space {R^n}. An /{F}-spline over
  {G}/ is a function {s} from {D} to {R^n}, such that {X*s = X*f_X}
  for every element {X} of {G}; where each {f_X} is some function of
  family {F}, and {X*f} denotes the restriction of {f} to {X}. Each
  partial function {X*f_X} is called a /piece/ of the spline.
  
  In particular, the restriction {D*f} of any function {f} of {F} to
  {D} is a /trivial/ {F}-spline over {G} (where all pieces are taken
  from the same function {f}).
  
  MULTI-DIMENSIONAL SPLINES
  
  In what follows, we consider the special case of meshes whose domain
  {D} is a subset of {R^d}, for some positive integer {d}, with the
  natural topology of the latter.
  
  For any {d} and any {i} in {0..d-1}, let {e^d_i} denote the unit
  vector of coordinate axis number {i} in {R^d}. */
  
typedef box_axis_t psp_axis_t; 
  /* Identifies a coordinate axis of a multidimensional 
    space, from 0 to {d-1} where {d} is the space's dimension. */

typedef box_axis_set_t psp_axis_set_t; 
  /* A packed subset of coordinate axes (each a {psp_axis_t}). */

typedef box_axis_index_t psp_axis_index_t; 
  /* Identifies a coordinate axis among a set of axes: 0 is the lowest. */

#define psp_NO_AXIS (BOX_NO_AXIS)
  /* A {psp_axis_t} value that means "no axis". */

#define psp_AXIS_NOT_FOUND (BOX_AXIS_NOT_FOUND)
  /* A {psp_axis_t} value that means "no such axis" */
  
psp_axis_set_t psp_axis_set_complement(psp_dim_t d, psp_axis_set_t A);
  /* The complement of {A} relative to the set {0..d-1}. */

/*
  CONTINUITY ORDER
    
  We say that a real-valued function {f} is /continuous to order {c}/,
  or /{c}-continuous/, with respect to a first-order differential
  operator {O}, if
  
    (1) {c} is {-1} and {f} is piecewise-continuous; or 
    
    (2) {c} is zero and {f} is continuous; or 

    (3) {c} is positive, {f} is continuous, and {O(f)} exists and is {(c-1)}-continuous.
        
  Note that a function with continuity {c} also has continuity {c'}
  for any {c'} in {-1..c}. */
  
typedef int8_t psp_cont_t;
  /* Continuity order of a function. The value {-1} means discontinuous.. */
 
/*
  MULTI-DIMENSIONAL CONTINUITY
  
  Let {c[0..d-1]} be a vector of {d} integers.  We say that a 
  function {f} is /continuous to order {c}/, or /{c}-continuous/,
  if {f} is continuosu to order {c[i]} with respect to 
  differentiation along each coordinate {i}.
  
  ----------------------------------------------------------------------
  
  SPLINE SPACES
  
  For any function family {F}, we denote by {\S(F,G)} the set of all
  {F}-splines over a mesh {G}. Note that if {F} is a
  finite-dimensional vector space, so is {\S(F,G)}.
  
  An /{F}-spline space over {G}/ is any vector space contained in {\S(F,G)}.

  CONTINUOUS SPLINE SUBSPACES
    
  We denote by {S_c(F,G)} the space of all splines in {S(F,G)} that
  have continuity order {c}. In particular, {S_0(F,G)} is the space of
  all continuous splines in {S(F,G)}.
  
  ----------------------------------------------------------------------
  
  MULTI-DIMENSIONAL POLYNOMIALS
  
  Let {g[0..d-1]} be a vector of {d} natural numbers. A /polynomial of
  degree {g}/ is a function {f} from {R^d} to {R} such that, for any {i} in
  {0..d-1} and any {x} in {R^d}, the function {t --> f(x + t*e^d_i)} 
  is a polynomial of degree {g[i]} on the parameter {t}.
  
  We denote by {P^g} the vector space of all such polynomials.
  POLYNOMIAL SPLINES
  
  A /polynomial spline of degree {g}/ is a spline on the family {P^g},
  i.e. an element of {S(P^g,G)} for some mesh {G}. */
  
typedef int8_t psp_degree_t;
  /* Degree of a polynomial, polynomial spline, etc. along some axis. */

/*
  
  CONTINUOUS POLYNOMIAL SPLINE SPACES
  
  We are particularly interested in the space {S_c(P^g,G)} 
  for some given mesh {G}; namely, the polynomial splines on {G}
  with pieces of degree {g[0..d-1]} with continuity {c[0..d-1]}. */

/*  
  SPLINE BASES

  An /{F}-spline basis over {G}/ is a basis for some {F}-spline space
  over {G}, that is, a set of linearly independent {F}-splines over
  {G}.
  
  If {J} is a basis for {F}, then the set of single-tile  
  
    {\set{ X*\phi : phi \in J \and X \in G }} 
    
  contains a basis for {\S(F,G)}. However, this may not be true for
  proper subspaces of {\S(F,G)}.
  
  ----------------------------------------------------------------------
   
  SPLINE SUPPORT
  
  The /{G}-support/ of a spline is the set of elements of {G} whose
  pieces not identically zero. Note that this set is not the same
  thing as the set of /points/ where the spline is nonzero.
 
  FINITE ELEMENTS

  A /finite element/ is a spline with finite {G}-support. (If {G} has
  finitely many elements, then technically any spline on {G} is a
  finite element; but we generally reserve the name for splines whose
  support is a relatively small subset of {G}.)

  ----------------------------------------------------------------------
   
  UNIDIMENSIONAL REGULAR SPLINES
  
  A /unidimensional regular spline/ is a polynomial spline defined on some
  unidimensional regular grid {G}, infinite or with circular topology.
  
  We will denote by {\S_c(P^g,G)} the space of all unidimensional
  splines over the grid {G}, whose pieces are polynomials of maximum
  degree {g} ({\geq 0}) and which are continuous to order {c} ({\geq
  -1}) everywhere.
  
  The {k}th derivative of a spline of order {c} is defined everywhere
  if {k <= c}; for {k > c}, it is defined only inside the cells, not
  at the grid vertices. In particular, a spline with continuity {c =
  -1} is undefined at the grid vertices. */

#endif
