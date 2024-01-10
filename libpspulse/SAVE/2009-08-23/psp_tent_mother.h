/* Mother tents for multidimensional toroidal regular grids. */
/* Last edited on 2009-08-23 19:37:59 by stolfi */

#ifndef psp_tent_H
#define psp_tent_H

#include <bz_basic.h>
#include <psp_pulse.h>
#include <psp_basic.h>
#include <bz_patch.h>
#include <vec.h>
#include <bool.h>

/*
  MOTHER TENTS
  
  A /mother tent/ is the product of {d} mother
  pulses {mp[0..d-1]}, each depending on only one coordinate: 
  
    { mt(x) = PROD { mp[i](x[i]) : i \in 0..d-1 } }
  
  Such a mother tent can be identified within a standard tent family
  {fam[0..d-1]} by the list {pix[0..d-1]}, where {pix[i]} is the index
  of the mother pulse {mp[i]} in the pulse family {fam[i]},
  
  The number of distinct mother tents defined by a tent family 
  {fam[0..d-1]} is, therefore,
  
    { nmt(fam) = PROD { nmp(fam[i]) : i \in 0..d-1 } }
    
  where {nmp(fam[i])} is the number of mother pulses in pulse
  family {fam[i]}. */
 
typedef uint32_t psp_tent_mother_count_t;
  /* A count of mother tents */

psp_tent_mother_count_t psp_tent_mother_count(psp_dim_t d, psp_pulse_family_t fam[]);
  /* The number of mother tents in the tent family {fam[]}. */

/*
  MOTHER TENT INDEX
  
  A /mother tent index/ {tix} is a single integer that conveniently
  combines the {d} mother pulse indices {pix[0..d-1]}, and therefore
  identifies a unique mother tent within the family.  The tent index
  is defined by the formula
  
    {tix = SUM { pix[i]*PROD {nmp(fam[j]) : 0 <= j < i} : i = 0..d-1 }
    
  where {nmp(fam[j])} is the number of mother pulses in the pulse 
  family {fam[j]}. */
  
typedef uint32_t psp_tent_mother_index_t;
  /* A packed vector of {psp_pulse_mother_index_t}s. */

#define psp_tent_mother_index_NONE (4294967295U)
  /* A {psp_tent_mother_index_t} value meaning `no mother tent index'. */

psp_tent_mother_index_t psp_tent_mother_index_max
  ( psp_dim_t d, 
    psp_pulse_family_t fam[]
  );
  /* The maximum packed index for any mother tent of family {fam[0..d-1]}. */

psp_tent_mother_index_t psp_tent_mother_index_pack
  ( psp_dim_t d,
    psp_pulse_family_t fam[],
    psp_pulse_mother_index_t pix[]
  );
  /* Packs pulse mother indices {pix[0..d-1]} into a tent mother index.  */

void psp_tent_mother_index_unpack
  ( psp_dim_t d, 
    psp_pulse_family_t fam[],
    psp_tent_mother_index_t tix,
    psp_pulse_mother_index_t pix[]
  );
  /* Unpacks a tent mother index {tix} into pulse mother indices {pix[0..d-1]}.  */

/*
  MOTHER TENT SUPPORT
  
  The support of a standard dyadic mother tent of family {fam[0..d-1]}
  and mother pulse indices {pix[0..d-1]} is a rectangular array of
  cells of the unit regular grid of {R^d}, which is the "Cartesian
  product" of the supports of its factor pulses {p[0..d-1]}.  */

void psp_tent_mother_supp_count
  ( psp_dim_t d,
    psp_pulse_family_t fam[],
    psp_pulse_mother_index_t pix[], 
    psp_grid_size_t msz[]
  );
  /* Saves in {msz[i]}, for {i} in {0..d-1}, the number of cells along
    axis in the support of the mother tent whose family is
    {fam[0..d-1]} and whose factor pulses have indices {pix[0..d-1]}.
    Namely sets {msz[i] = psp_pulse_mother_supp_count(fam[i],pix[i])}, 
    for each {i}. */ 

#endif
