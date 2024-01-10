/* Finite elements for spline spaces on multidimensional toroidal uniform grids. */
/* Last edited on 2009-08-23 19:57:16 by stolfi */

#ifndef psp_tent_H
#define psp_tent_H

#include <bz_basic.h>
#include <psp_pulse.h>
#include <psp_basic.h>
#include <bz_patch.h>
#include <vec.h>
#include <bool.h>

/*
  TENTS
  
  A /tent/ is the multi-dimensional analog of a pulse: it is a finite
  element (a polynomial spline with finite support) defined on a
  multidimensional regular circular grid.
  
  STANDARD TENT FAMILIES
  
  Here we define a number of /standard tent families/ that are used as
  bases of spline spaces defined on regular toroidal multidimensional
  grids. For the time being, a family of standard tents is defined by
  the dimension {d} of the domain and {d} standard pulse families
  {fam[0..d-1]}. Each member of that tent family is the product of {d}
  pulses, each taken from a standard pulse family and depending on
  only one coordinate of the domain. More precisely, the value of a
  standard {d}-dimensional tent {t} at a generic point {x} of {R^d} is
  
    { t(x) = PROD { p[i](x[i]) : i \in 0..d-1 } }
  
  where {p[0..d-1]} are standard pulses.
  
  In what follows, /tent/ will usually mean a standard tent. Also
  {fam[i] = (pkind[i],c[i],g[i])} is the family of the factor pulse
  {p[i]}; {pos[i] = p[i].pos} is the displacement of that pulse along
  axis {i}; {pix[i] = p[i].pix} is its mother pulse index; and {gsz[i]
  = p[i].gsz} is the size (number of cells) of the grid along axis {i}.
  
  STANDARD TENT EVALUATION
  
  The following procedures compute the value and partial derivatives of
  a dyadic tent {t} at a point {x} of {R^d}. They properly handle
  wrap-around and self-overlapping of the mother pulses.
  
  The results are stored in {f[0..N-1]}, for some {N}: the function
  value is stored at {f[0]}, followed by the other derivatives, as
  requested.
  
  The root cell {R} defaults to {(0_1)^d} if {R == NULL}. */

void psp_tent_eval
  ( psp_dim_t d,
    psp_pulse_family_t fam[], 
    psp_pulse_t p[],
    interval_t R[],
    double x[],
    psp_degree_t ord[],
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

void psp_tent_eval_total
  ( psp_dim_t d,
    psp_pulse_family_t fam[], 
    psp_pulse_t p[],
    interval_t R[],
    double x[],
    psp_degree_t ord,
    double f[]
  );
  /* Same as {psp_tent_eval}, but computes only derivatives up to 
    *total* order {ord} in each coordinate {x[i]}. The derivatives are
    stored in {f[0..N-1]}, where {N = choose(ord+d,d) }, sorted
    by total order and then lexicographically. E.g., for {d=3,ord=2},
    the {f} vector will contain the {choose(5,3) = 10} derivatives

       { (w,wx,wy,wz,wxx,wxy,wyy,wxz,wyz,wzz) }. */
    
/*
  STANDARD TENT SUPPORT
  
  Let {t} be a standard dyadic tent of rank {r}, with tent family
  {fam[0..d-1]} and mother pulse indices {pix[0..d-1]}. The support of
  {t} is a cell pack of rank {r}, which is the "Cartesian product" of
  the supports of its factor pulses {p[0..d-1]}.
  
  For large enough grids, the number {psz[i]} of cells in {t}'s
  support along any axis {i} will be {msz[i] =
  psp_pulse_mother_supp_count(fam[i],pix[i])}. However, if the grid
  size {gsz[i]} along axis {i} is less than {msz[i]}, the toroidal
  topology will cause the pulse to wrap around the boundary of {D} and
  overlap itself.  Note that the wrap-around may happen only along some
  axes.  In general,
  
    {psz[i] = psp_pulse_supp_count(fam[i],pix[i],gsz[i])}
  
  ANCHOR CELL OF A TENT
  
  The /anchor cell/ of a tent {t} is the cell in {t}'s
  support that corresponds to the lowest cell of the 
  mother tent's support. The cells of pack {C} is linearized
  in such a way that its first cell is the anchor cell.
  
  The position of the anchor cell in level {r} is {pos[0..d-1]}, where
  {pos[i]} is the shift of pulse {p[i]} along axis {i}.
    
  The tent's support {B} then spans the the interval
  {pos[i]..pos[i]+psz[i]-1} along each axis {i}, modulo {gsz[i]}.
*/

/* PRINTOUT */
  
void psp_tent_print
  ( FILE *wr, 
    psp_dim_t d,
    psp_pulse_family_t fam[], 
    psp_pulse_t p[]
  );
  /* Prints the tent whose factors are pulses {p[0..d-1]} from
    families {fam[0..d-1]}, in the format "{pulse[0]}*...{pulse[d-1]}"
    where {pulse[i]} is the output of {psp_pulse_print}. */

#endif
