/* H-pulses - Hermite-like univariate finite elements. */
/* Last edited on 2011-09-18 13:50:10 by stolfilocal */

#ifndef psp_pulse_mother_H_H
#define psp_pulse_mother_H_H

#include <psp_pulse.h>

/* 
  H-PULSES
  
  The dyadic H-pulses are defined only for {c\geq 0}, {g = 2*c+1}.
  
  There are exactly {c+1} mother H-pulses {w^g_c[pix]}. Their support is
  the interval {(0 _ 2)}. The mother pulse of index {pix} has all
  derivatives up to order {c} equal to zero at integer {z} values,
  except for the derivative of order {pix} at {z = 1}, which is equal to
  {choose(g,pix)}.
  
  The mother pulse is symmetric or anti-symmetric around {z = 1}, 
  depending on the  parity of {pix}; so, for {z} in {(0,1)}, 

    { w^g_c[pix](z) = (-1)^pix w^g_c[pix](2-z) }
  
  */

unsigned int psp_pulse_mother_H_count (psp_cont_t c, psp_degree_t g);
  /* Returns the number of mother H-pulses with continuity {c}
    and degree {g}; namely, {c + 1} if {c >= 0} and 
    {g = 2*c + 1}, and 0 otherwise. */

psp_grid_size_t psp_pulse_mother_H_supp_count
  ( psp_cont_t c,
    psp_degree_t g,
    psp_pulse_mother_index_t pix
  ); 
  /* Returns the number of unit-grid cells that comprise the support
    of the mother H-pulse with given {c,g,pix}. Fails if such H-pulses
    aren't defined or aren't supported by this module. */

void psp_pulse_mother_H_to_bezier
  ( psp_cont_t c, 
    psp_degree_t g,
    psp_pulse_mother_index_t pix,
    psp_grid_pos_t x,
    double *bz
  );
  /* Stores into {bz[0..2*c+1]} the Bézier coefficients of the mother
    H-pulse with continuity {c}, degree {g}, and mother index {pix}, for
    the interval {(x _ x+1)}. Fails if such pulse is not defined or
    supported by this module. The coefficients are nonzero only if {x}
    is 0 or 1. */

#endif
