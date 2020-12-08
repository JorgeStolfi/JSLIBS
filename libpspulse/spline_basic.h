#ifndef spline_basic_H
#define spline_basic_H

/* Basic definitions for splines on generic meshes. */
/* Last edited on 2011-09-19 00:57:07 by stolfilocal */

#define _GNU_SOURCE
#include <stdint.h>

/*
  GENERAL MESHES
  
  A /mesh/ is a partition {G} of some topological space {D} (the /domain/ of
  {G}) into a finite set of subsets (the /items/ of {G}) such that (1) each item
  is a topological open {k}-ball of {D}, for some finite {k}; and (2) the
  topological closure in {D} of any item of {G} is the union of finitely many
  items of {G}.
  
  The items with maximum dimension are the /cells/ of the mesh.
  Items with dimension {k = 0,1,2,3} are called /vertices/, /edges/, 
  /panes/, and /blocks/, respectively.
  
  A /submesh/ is any subset of a mesh; the /domain/ of the submesh is 
  is the union of its elements.
  
  MESH REFINEMENT
  
  Amesg {H} is said to be a /refinement/ of a mesh {G} if every item
  of {G} is the union of fintely many items of {H}.
  
  ----------------------------------------------------------------------
  
  GENERAL SPLINES
  
  If {f} is any function and {X} is a subset of its domain, we will denote by
  {X*f} the restriction of {f} to {X}.
  
  Let {G} be a mesh with domain {D}, and {F} be a family of functions from some
  superset of {D} to some vector space {V}. An /{F}-spline over {G}/ is a
  function {s} from {D} to {V}, such that {X*s = X*f_X} for every item {X} of
  {G}; where each {f_X} is some function of family {F}. Each partial function
  {X*f_X} is called a /piece/ of the spline.
  
  If all pieces of {f} are restrictions of the same function of {F},
  then {f} is said to be /trivial/.
  
  Note that any {F}-spline over any mesh {G} is also an {F}-spline
  over any refinement of {G}.  
  
  Let {f} be an {F}-spline over a mesh {G}. Let {H} be a submesh of {G}, and {E} the domain
  of {H}. Note that the restriction {E*f} of {f} to {E} is an {F}-spline over {H}.
    
  ----------------------------------------------------------------------
  
  SPLINE SPACES
  
  We denote by {\S(F,G)} the set of all {F}-splines over {G}. Note that if {F}
  is a finite-dimensional vector space, so is {\S(F,G)}.
  
  An /{F}-spline space over {G}/ is any vector space contained in {\S(F,G)}.
 
  SPLINE BASES

  An /{F}-spline basis over {G}/ is a basis for some {F}-spline space
  over {G}, that is, a set of linearly independent {F}-splines over
  {G}.
  
  If {J} is a basis for {F}, then the set of single-tile splines 
  
    {\set{ X*\phi : phi \in J \and X \in G }} 
    
  contains a basis for {\S(F,G)}. However, this may not be true for
  proper subspaces of {\S(F,G)}.
   
  ----------------------------------------------------------------------
  
  SPLINE SUPPORT
  
  The /{G}-support/ of a spline is the set of items of its grid {G} whose
  pieces not identically zero. Note that this set is not the same thing as the
  set of /points/ where the spline is nonzero.
  
  FINITE ELEMENTS

  A /finite element/ is a spline with finite {G}-support. (If {G} has
  finitely many elements, then technically any spline on {G} is a
  finite element; but we generally reserve the name for splines whose
  support is a relatively small subset of {G}.)
  
  FINITE ELEMENT BASES
  
  A /finite element basis/ for any spline space {S} on a grid {G} is a 
  basis for {S} consisting of finite elements.
  
  ----------------------------------------------------------------------
  
  MONOMIALS
  
  Let {D} be a {d}-dimensional subset of some Cartesian space {\RR^d}. A
  /monomial function/ (or just /monomial/) /on {D}/ is a function {f} from {D}
  to {\RR}, whose value {f(x)} at any point {x} of {D} is given by a product of
  non-negative integer powers of the coordinates,
  
    {f(x) = x[0]^{e[0]} * x[1]^{e[1]} * ... * x[d-1]^{e[d-1]}}
    
  The vector {e[0..d-1]} is the /degree vector/ of the monomial {f}. The /total
  degree/ of {f} is the sum of all {e[i]}. The total degree is zero if and only if
  all {e[i]} are zero, in which case {f} is the constant function equal to 1
  everywhere.
  POLYNOMIALS
  
  A /polynomial function on {D}/ is a function {f} from {D} to some linear space
  {V}, which is a linear combination of monomials with fixed coefficients from
  {V}.
  
  The /total degree/ of {f} is the largest total degree {t} among its 
  monomials with nonzero coefficients.  Its /degree vector/ 
  is the vector {m[0..d-1]} such that {m[i]} is the maximum 
  exponent of coordinate {x[i]} among those monomials.  If {f}
  is identically zero, {t} and all {m[i]} are {-oo} by convention.
  
  We will denote by {\P(D,V)} the set of all polynomial functions from {D} to {V}. For
  any natural number {g}, we will denote by {\P^g(D,V)} the subset of all
  polynomials in {\P(D,V)} with total degree {t} at most {g}. For any vector
  {g[0..d-1]} of natural numbers, we will denote by {\P^g(D,V)} the subset of all
  polynomials in {\P(D,V)} such that their degree vector {m[0..d-1]} satisfies
  {m[i] <= g[i]} for all {i}.
  
  In all these spaces we will drop {D} when clear from the context, and {V} when 
  it is the reals {\RR}. */
  
typedef int8_t spline_degree_t;
  /* Degree of a monomial, polynomial, spline, etc.. */

/* 
  CONTINUITY ORDER
  
  A function {f} is said to be /piecewise continuous/ on a mesh {G} if, for each
  item {X} of {G}, it is continuous within {X}, and tends to a finite limit
  whenever the argument tends to a point on the boundary of {X} by any path
  within {X}.
  
  Let {f} be a function defined on a {d}-dimensional subset of {\RR^d}, and {G}
  a mesh on {D}. We say that {f} is /continuous to order {c} on {G}/ if
  
    (1) {c = -1}, and {f} is piecewise continuous on {G};  or
    
    (2) {c = 0}, and {f} is continuous on {D}, differentiable within each item
      of {G}, and its gradient is piecewise continuous on {G}; or

    (2) {c > 0}, and {f} is continuous and differentiable at any point {p} of
      {D}, and its gradient is continuous to order {c-1} on {G}.
      
   Note that a function with continuity {c} also has continuity {c'} 
   for any {c'} in {-1..c}.
  */

typedef int8_t spline_cont_t;
  /* Continuity order of a function {f}. */

#endif
