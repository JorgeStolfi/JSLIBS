#ifndef dg_spline_H
#define dg_spline_H

/* Polynomial splines in arbitrary dimension */
/* Last edited on 2007-10-12 20:06:08 by stolfi */

#include <stdint.h>

/*
  GENERAL SPLINES
  
  A /{d}-dimensional mesh/ is a collection of pairwise disjoint
  subsets of {R^d}, the /tiles/ of the mesh.
  
  Let {L} be a {d}-dimensional mesh, and {F} be a family of functions
  from {R^d} to {R}. An /{F}-spline over {L}/ is a function {s} from
  {D = \U L} to {R}, such that {X*s = X*f_X} for every tile {X} of {L};
  where each {f_X} is some function of {F}, and {X*f} denotes the
  restriction of {f} to {X}. Each partial function {X*f_X} is called a
  /piece/ of the spline.
  
  In particular, the restriction {D*f} of any function {f} of {F} to
  {D = \U L} is a /trivial/ {F}-spline over {L} (where all pieces use
  the same function {f}).
  
  DYADIC SPLINES
  
  A spline is /dyadic/ if its mesh {L} is the set of leaf faces of
  some finite dyadic grid {G}.
  
  SPLINE SUPPORT
  
  The /{L}-support/ of a spline is the set of tiles of {L} whose
  pieces are not identically zero.  Note that this set is generally
  larger than the set of /points/ where the spline is nonzero.
  
  FINITE ELEMENTS

  A /finite element/ is a spline with finite, and generally small, {L}-support.
  
  SPLINE SPACES
  
  For any function family {F}, we denote by {\S(F,L)} the set of all
  {F}-splines over a mesh {L}. Note that if {F} is a
  finite-dimensional vector space, so is {\S(F,L)}.
  
  An /{F}-spline space over {L}/ is any vector space contained in {\S(F,L)}.
  
  SPLINE BASES

  An /{F}-spline basis over {L}/ is a basis for some {F}-spline space over {L},
  that is, a set of linearly independent {F}-splines over {L}.
  
  If {J} is a basis for {F}, then the set of single-tile  
  
    {\set{ X*\phi : phi \in J \and X \in L }} 
    
  contains a basis for {\S(F,L)}. However, this set may not be true 
  for proper subspaces of {\S(F,L)}.
  
  CONTINUITY ORDERS
  
  A real function {f} defined on some subset {X} of {R} is said to
  /have continuity order {c}/ (or just /continuity {c}/) iff
  
    (0) {c == -1} and {f} is piecewise-continuous in {X}, or
    
    (1) {c == 0} and {f} is continuous in {X}, or
    
    (2) {c > 0}, {f} is differentiable in {X}, and its derivative
        has continuity {c-1} in {X}. 
        
  Note that a function with continuity {c} also has continuity {c'}
  for any {c'} in  {-1..c}. */
  
typedef int8_t dg_cont_t;
  /* Continuity order of a function. The value {-1} means discontinuous.. */
 
/*
  CONTINUOUS SPLINE SPACES
    
  Let {e^d_i} denote the unit vector of {R^d} along axis {i}; that is,
  a vector such that {e^d_i[j]} is 1 if {i==j}, and 0 otherwise.
  
  Let {c[0..d-1]} be a vector of continuity orders. We say that a
  function {f} from some subset {D} of {R^d} to {R} has /continuity {c}/ if, for each
  axis {i} in {0..d-1}, and any {x} in {R^d}, the function {t --> f(x
  + t*e^d_i)} has continuity {c[i]}.
  
  We denote by {S_c(F,L)} the subspace of {S(F,L)} consisting of all
  {F}-splines over {L} that have continuity {c}.
  
  MULTIVARIATE POLYNOMIALS
  
  Let {g[0..d-1]} be a vector of {d} natural numbers. A /polynomial of
  degree {g}/ is a function {f} from {R^d} to {R} such that, for any {i} in
  {0..d-1} and any {x} in {R^d}, the function {t --> f(x + t*e^d_i)} 
  is a polynomial of degree {g[i]} on the parameter {t}.
  
  We denote by {P^g} the vector space of all such polynomials.
  
  POLYNOMIAL SPLINES
  
  A /polynomial spline of degree {g}/ is a spline on the family {P^g},
  i.e. an element of {S(P^g,L)} for some mesh {L}.
  
  Splines can be assumed to be polynomial unless said otherwise. 
  
  CONTINUOUS POLYNOMIAL SPLINE SPACES
  
  We are particularly interested in the space {S_c(P^g,G)} where {G} is a 
  dyadic grid; namely, the /dyadic polynomial splines with continuity {c[0..d-1]}
  and degree {g[0..d-1]}/ */
  
typedef int8_t dg_degree_t;
  /* Degree of a polynomial or spline, order of a derivative, etc.. */

#endif
