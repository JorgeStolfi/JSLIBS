/* Tent bases for polynomial splines defined on irregular dyadic grids */
/* Last edited on 2009-05-18 14:10:01 by stolfi */

#ifndef dg_tent_basis_H
#define dg_tent_basis_H

#include <bz_basic.h>
#include <dg_tent.h>
#include <dg_tree.h>
#include <udg_pulse.h>
#include <vec.h>

/*
  TENT BASES FOR CONTINUOUS SPLINES
  
  Let {fam[0..d-1]} be an arbitrary vector of pulse families, and let
  {fam[i] = (pkind[i],c[i],g[i])}. If {G} is a finite dyadic grid, the
  space {S_c(P^g,G)} admits /a tent basis of family {fam[]}/,
  a basis where all elements are tents of that family.
  
  FROM TEMPLATES TO BASES
  
  To build a tent basis of family {fam[0..d-1]} for the space
  {S_c(P^g,G)}, we must find a complete and linearly independent set
  of tents of that family whose supports consist of cells of {G}.
  
  Thus, for each mother tent index vector {pix[0..d-1]} that is valid
  for that family we must enumerate all cell packs that fit the
  template {msz[0..d-1]}, where
  
    {msz[i] = udg_pulse_mother_supp_size(fam[i],pix[i])}
    
  For every such cell pack {C}, we then take the dyadic tent defined
  by said mother tent, scaled and shifted so that its support is {C}.

  The result is a set of tents that generate the space {S_c(P^g,G)}.
  Howeer, that set usually contains linearly dependent subsets,
  which should be eliminated in order to create a basis.  This 
  module allows the client to choose between two linearly independent
  subsets: the /inferior basis/ and the /superior basis/.
  
  INFERIOR BASIS
  
  To build the inferior basis, one uses an instance {C} of the 
  template if and only if it contains at least one leaf cell.
  
  SUPERIOR BASIS
  
  To build the superior basis, one uses an instance {C} of the 
  template if and only if its lowest cell {C[0]} has even index. */
  
dg_tent_vec_t dg_tent_basis_get
  ( mdg_dim_t d, 
    udg_pulse_family_t fam[], 
    dg_tree_t G,
    bool_t superior
  );
  /* Returns a set of {d}-dimensional tents, consisting of all
    possible products of pulses from family {fam[i]} along each axis
    {i} in {0..d}.
    
    If {superior} is TRUE, returns the `superior' basis, else
    returns the `inferior' basis. */

/*
  COMPLETE VERTICES

  In particular, a rank-{r} cell {C} of {G} fits the template
  {(2,2,..2)} if and only if every rank-{r} dyadic cell that is
  incident to its upper corner {v} is a cell of {G}. Thus, for mother
  tents whose support size is that template, we need to locate the
  /complete vertices/ of {G} -- namely those vertices which are
  incident to {2^d} cells of the same level (counting repetitions). */

#endif
