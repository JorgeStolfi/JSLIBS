/* Finite elements for multidimensional dyadic spline spaces. */
/* Last edited on 2009-05-18 14:09:08 by stolfi */

#ifndef dg_tent_H
#define dg_tent_H

#include <mdg_grid.h>
#include <bz_basic.h>
#include <udg_pulse.h>
#include <dg_spline.h>
#include <bz_patch.h>
#include <vec.h>
#include <bool.h>

/*
  DYADIC TENTS
  
  A /dyadic tent/ is a polynomial spline (typically with small and compact support)
  defined on some level {r} of the {d}-dimensional dyadic multigrid {G*}.
  
  STANDARD DYADIC TENTS
  
  For multidimensional grids, the finite-element bases that we
  consider here consist of /standard dyadic finite tents/. Indeed, a
  standard dyadic tent is the product of {d} standard dyadic pulses,
  {b[0..d-1]}, each depending on only one coordinate:
  
    { t(x) = PROD { b[i](x[i]) : i \in 0..d-1 } }
  
  In what follows, /dyadic tent/ will usually mean a standard dyadic
  tent. Also {fam[i] = (pkind[i],c[i],g[i])} is the family of the
  factor pulse {b[i]}; {pos[i]} is the displacement of that pulse;
  {pix[i]} is its mother pulse index; and {gsz[i]} is the size along
  axis {i} of the grid where the tent lives.
  
  While the parameters {pkind[i],c[i],g[i],pix[i],pos[i]} can be chosen
  independently for each {i}, the grid size of {b[i]} is determined by
  the level {r} of the grid, i.e. {gsz[i] = mdg_grid_size(d,r,i)}.
  
  STANDARD TENT EVALUATION
  
  The following procedures compute the value and partial derivatives of
  a dyadic tent {t} at a point {x} of {R^d}. They properly handle
  wrap-around and self-overlapping of the mother pulses.
  
  The results are stored in {f[0..N-1]}, for some {N}: the function
  value is stored at {f[0]}, followed by the other derivatives, as
  requested.
  
  The root cell {R} defaults to {(0_1)^d} if {R == NULL}. */

void dg_tent_eval
  ( mdg_dim_t d,
    udg_pulse_family_t fam[], 
    udg_pulse_t p[],
    interval_t R[],
    double x[],
    dg_degree_t ord[],
    double f[]
  );
  /* Evaluates a tent function and all its partial derivatives up to 
    order {ord} at a point {x[0..d-1]}.
    
    Computes all derivatives up to order {ord[i]} in each
    coordinate {x[i]}. The derivatives are stored in {f[0..N-1]},
    where  {N = PROD { ord[i]+1 : i = 0..d-1 } },
    in ``column by column'' order. E.g., for {d=3,ord=(2,1,2)},
    the {f} vector will contain the {3*2*3 = 18} derivatives
        
      { (w,wx,wxx,wy,wxy,wxxy,wz,wxz,wxxz,wyz,wxyz.wxxyz,wzz, ... wxxyzz) }
      
    where {wxxy} means {w} differentiated twice on {x[1]} and once on {x[2]}. */

void dg_tent_eval_total
  ( mdg_dim_t d,
    udg_pulse_family_t fam[], 
    udg_pulse_t p[],
    interval_t R[],
    double x[],
    dg_degree_t ord,
    double f[]
  );
  /* Same as {dg_tent_eval}, but computes only derivatives up to 
    *total* order {ord} in each coordinate {x[i]}. The derivatives are
    stored in {f[0..N-1]}, where {N = choose(ord+d,d) }, sorted
    by total order and then lexicographically. E.g., for {d=3,ord=2},
    the {f} vector will contain the {choose(5,3) = 10} derivatives

       { (w,wx,wy,wz,wxx,wxy,wyy,wxz,wyz,wzz) }. */
  
/*
  STANDARD TENT FAMILIES
  
  A /standard tent family/ is a vector {fam[0..d-1]} of pulse
  families.
  
  MOTHER TENTS
  
  A /mother tent/ is a spline {mt(x)} defined on the infinite
  {d}-dimensional unit regular grid as the product of {d} mother
  pulses {mb[0..d-1]}, each depending on only one coordinate: 
  
    { mt(x) = PROD { mb[i](x[i]) : i \in 0..d-1 } }
  
  Such a mother tent can be identified within a standard tent family
  {fam[0..d-1]} by the list {pix[0..d-1]}, where {pix[i]} is the index
  of the mother pulse {mb[i]} in the pulse family {fam[i]},
  
  The number of distinct mother tents defined by a tent family 
  {fam[0..d-1]} is, therefore,
  
    { nmt(fam) = PROD { nmp(fam[i]) : i \in 0..d-1 } }
    
  where {nmp(fam[i])} is the number of mother pulses in pulse
  family {fam[i]}. */
 
typedef uint32_t dg_tent_mother_count_t;
  /* A count of mother tents */

dg_tent_mother_count_t dg_tent_mother_count(mdg_dim_t d, udg_pulse_family_t fam[]);
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
  
typedef uint32_t dg_tent_mother_index_t;
  /* A packed vector of {udg_pulse_mother_index_t}s. */

#define dg_tent_mother_index_NONE (4294967295U)
  /* A {dg_tent_mother_index_t} value meaning `no mother tent index'. */

dg_tent_mother_index_t dg_tent_mother_index_max(mdg_dim_t d, udg_pulse_family_t fam[]);
  /* The maximum packed index for any mother tent of family {fam[0..d-1]}. */

dg_tent_mother_index_t dg_tent_mother_index_pack
  ( mdg_dim_t d,
    udg_pulse_family_t fam[],
    udg_pulse_mother_index_t pix[]
  );
  /* Packs pulse mother indices {pix[0..d-1]} into a tent mother index.  */

void dg_tent_mother_index_unpack
  ( mdg_dim_t d, 
    udg_pulse_family_t fam[],
    dg_tent_mother_index_t tix,
    udg_pulse_mother_index_t pix[]
  );
  /* Unpacks a tent mother index {tix} into pulse mother indices {pix[0..d-1]}.  */

/*
  MOTHER TENT SUPPORT
  
  The support of a standard dyadic mother tent of family {fam[0..d-1]}
  and mother pulse indices {pix[0..d-1]} is a rectangular array of
  cells of the unit regular grid of {R^d}, which is the "Cartesian
  product" of the supports of its factor pulses {b[0..d-1]}.  */

void dg_tent_mother_supp_count
  ( mdg_dim_t d,
    udg_pulse_family_t fam[],
    udg_pulse_mother_index_t pix[], 
    mdg_grid_size_t msz[]
  );
  /* Saves in {msz[i]}, for {i} in {0..d-1}, the number of cells along
    axis in the support of the mother tent whose family is
    {fam[0..d-1]} and whose factor pulses have indices {pix[0..d-1]}.
    Namely sets {msz[i] = udg_pulse_mother_supp_count(fam[i],pix[i])}, 
    for each {i}. */ 
    
/*
  STANDARD TENT SUPPORT
  
  Let {t} be a standard dyadic tent of rank {r}, with tent family
  {fam[0..d-1]} and mother pulse indices {pix[0..d-1]}. The support of
  {t} is a cell pack of rank {r}, which is the "Cartesian product" of
  the supports of its factor pulses {b[0..d-1]}.
  
  For large enough grids, the number {psz[i]} of cells in {t}'s
  support along any axis {i} will be {msz[i] =
  udg_pulse_mother_supp_count(fam[i],pix[i])}. However, if the grid
  size {gsz[i]} along axis {i} is less than {msz[i]}, the toroidal
  topology will cause the pulse to wrap around the boundary of {D} and
  overlap itself.  Note that the wrap-around may happen only along some
  axes.  In general,
  
    {psz[i] = udg_pulse_supp_count(fam[i],pix[i],gsz[i])}
  
  ANCHOR CELL OF A TENT
  
  The /anchor cell/ of a tent {t} is the cell in {t}'s
  support that corresponds to the lowest cell of the 
  mother tent's support. The cells of pack {C} is linearized
  in such a way that its first cell is the anchor cell.
  
  The position of the anchor cell in level {r} is {pos[0..d-1]}, where
  {pos[i]} is the shift of pulse {b[i]} along axis {i}.
    
  The tent's support {B} then spans the the interval
  {pos[i]..pos[i]+psz[i]-1} along each axis {i}, modulo {gsz[i]}.

  TENT REPRESENTATION
  
  In conclusion, a specific tent within a given tent family can be specified
  completely by two integers: its mother tent index {tix}, and the cell
  index {ck} of its anchor cell.  Note that {ck} determines the level {r},
  and hence the grid sizes {gsz[0..d-1]}, as well as the shifts {pos[0..d-1]}. */

typedef struct dg_tent_t  /* Identifies a tent within a tent family. */
  { mdg_cell_index_t ck;         /* Index of anchor cell. */
    dg_tent_mother_index_t tix;  /* Packed mother pulse indices. */
  } dg_tent_t;
  /* Identifies (within some family fixed by context) the tent whose
    mother tent index is {tix} and whose anchor cell has cell index
    {ck}.  */

dg_tent_t dg_tent_pack
  ( mdg_dim_t d, 
    udg_pulse_family_t fam[],
    udg_pulse_t p[]
  );
  /* Returns the packed description of the tent formed by pulses
    {p[0..d-1]} from families {fam[0..d-1]}, respectively. */

void dg_tent_unpack
  ( mdg_dim_t d, 
    udg_pulse_family_t fam[],
    dg_tent_t *t,
    udg_pulse_t p[]
  );
  /* Unpacks a packed tent description {t} into its factor pulses
    {p[0..d-1]}, assuming they belong to families {fam[0..d-1]}. */

/* 
  ENUMERATING TENTS */

vec_typedef(dg_tent_vec_t,dg_tent_vec,dg_tent_t);
  /* A list of {dg_tent}s. */

dg_tent_vec_t dg_tents_from_cells
  ( mdg_dim_t d, 
    dg_tent_mother_index_t tix,
    mdg_cell_index_vec_t ck
  );
  /* Returns a list of {dg_tent_t} {t} obtained by combining tent
    index {tix} with all cell indices {ck[0..]} */

/* PRINTOUT: */
  
void dg_tent_print
  ( FILE *wr, 
    mdg_dim_t d,
    udg_pulse_family_t fam[], 
    udg_pulse_t p[]
  );
  /* Prints the tent whose factors are pulses {p[0..d-1]} from
    families {fam[0..d-1]}, in the format "{pulse[0]}*...{pulse[d-1]}"
    where {pulse[i]} is the output of {udg_pulse_print}. */

void dg_tent_packed_print
  ( FILE *wr, 
    mdg_dim_t d,
    udg_pulse_family_t fam[], 
    dg_tent_t *t
  );
  /* Like {dg_tent_print} but extracts the tent's info
    from the packed tent {t}. */

#endif
