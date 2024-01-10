#ifndef psp_basic_H
#define psp_basic_H

/* Basic definitions for polynomial splines on generic meshes. */
/* Last edited on 2009-08-23 13:19:35 by stolfi */

#include <stdint.h>

/*
  GENERAL MESHES
  
  A /mesh/ is a finite partition {L} of some topological space {D} (the
  /domain/) such that (1) each element of {L} is a topological open
  {k}-ball of {D}, for some {k}; and (2) the topological closure in
  {D} of any element of {L} is the union of finitely many elements of
  {L}. 
  
  The elements with maximum dimension are the /cells/ of the mesh.
  Elements with dimension {k = 0,1,2,3} are called /vertices/, /edges/, 
  /panes/, and /blocks/, respectively.
  
  GENERAL SPLINES
  
  Let {L} be a mesh with domain {D}, and {F} be a family of functions
  from {D} to some vector space {R^n}. An /{F}-spline over {L}/ is a
  function {s} from {D} to {R^n}, such that {X*s = X*f_X} for every
  element {X} of {L}; where each {f_X} is some function of family {F},
  and {X*f} denotes the restriction of {f} to {X}. Each partial
  function {X*f_X} is called a /piece/ of the spline.
  
  In particular, the restriction {D*f} of any function {f} of {F} to
  {D} is a /trivial/ {F}-spline over {L} (where all pieces are taken
  from the same function {f}). */
  
typedef int8_t psp_ddim_t;
  /* Dimension of the domain of a spline. */

typedef int8_t psp_rdim_t; 
  /* Dimension of range of a spline. */

/* 
  SPLINE SPACES
  
  For any function family {F}, we denote by {\S(F,L)} the set of all
  {F}-splines over a mesh {L}. Note that if {F} is a
  finite-dimensional vector space, so is {\S(F,L)}.
  
  An /{F}-spline space over {L}/ is any vector space contained in {\S(F,L)}.

  ----------------------------------------------------------------------

  CONTINUITY ORDER
    
  We say that a real-valued function {f} is /continuous to order {c}/,
  or /{c}-continuous/, with respect to a first-order differential
  operator {O}, if
  
    (1) {c} is negative and {f} is piecewise-continuous; or 
    
    (2) {c} is zero and {f} is continuous; or 
    
    (3) {c} is positive, {f} is continuous, and {O(f)} exists and is {(c-1)}-continuous.
        
  Note that a function with continuity {c} also has continuity {c'}
  for any {c'} in {-1..c}. */
  
typedef int8_t psp_cont_t;
  /* Continuity order of a function. The value {-1} means discontinuous.. */
 
/*
  MULTI-DIMENSIONAL SPLINES
  
  In what follows, we consider the special case when {D} is a subset
  of {R^d}, for some positive integer {d}, with the natural topology
  of the latter.
  
  For any and any {i} in {0..d-1}, let {e^d_i} denote
  the unit vector of coordinate axis number {i} in {R^d}.
  
  MULTI-DIMENSIONAL CONTINUITY
  
  Let {c[0..d-1]} be a vector of {d} integers.  We say that a 
  function {f} is /continuous to order {c}/, or /{c}-continuous/,
  if {f} is continuosu to order {c[i]} with respect to 
  differentiation along coordinate {i}.
  
  CONTINUOUS SPLINE SUBSPACES
    
  We denote by {S_c(F,L)} the space of all splines in {S(F,L)} that
  have continuity order {c}. In particular, {S_0(F,L)} is the space of
  all continuous splines in {S(F,L)}.
  
  MULTI-DIMENSIONAL POLYNOMIALS
  
  Let {g[0..d-1]} be a vector of {d} natural numbers. A /polynomial of
  degree {g}/ is a function {f} from {R^d} to {R} such that, for any {i} in
  {0..d-1} and any {x} in {R^d}, the function {t --> f(x + t*e^d_i)} 
  is a polynomial of degree {g[i]} on the parameter {t}.
  
  We denote by {P^g} the vector space of all such polynomials.
  
  POLYNOMIAL SPLINES
  
  A /polynomial spline of degree {g}/ is a spline on the family {P^g},
  i.e. an element of {S(P^g,L)} for some mesh {L}.
  
  CONTINUOUS POLYNOMIAL SPLINE SPACES
  
  We are particularly interested in the space {S_c(P^g,L)} 
  for some given mesh {L}; namely, the polynomial splines on {L}
  with pieces of degree {g[0..d-1]} with continuity {c[0..d-1]}. */
  
typedef int8_t psp_degree_t;
  /* Degree of a polynomial or spline etc.. */

/*  
  SPLINE BASES

  An /{F}-spline basis over {L}/ is a basis for some {F}-spline space
  over {L}, that is, a set of linearly independent {F}-splines over
  {L}.
  
  If {J} is a basis for {F}, then the set of single-tile  
  
    {\set{ X*\phi : phi \in J \and X \in L }} 
    
  contains a basis for {\S(F,L)}. However, this may not be true for
  proper subspaces of {\S(F,L)}.
  
  SPLINE SUPPORT
  
  The /{L}-support/ of a spline is the set of tiles of {L} whose
  pieces are not identically zero.  Note that this set is not the same
  thing as the set of /points/ where the spline is nonzero.
  
  FINITE ELEMENTS

  A /finite element/ is a spline with finite, and generally small, {L}-support.  */


#endif
