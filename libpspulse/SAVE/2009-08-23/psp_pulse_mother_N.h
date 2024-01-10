/* N-pulses - Compact unit-partition finite elements. */
/* Last edited on 2009-08-23 18:25:19 by stolfi */

#ifndef psp_pulse_mother_N_H
#define psp_pulse_mother_N_H

#include <bz_basic.h>

#include <psp_pulse.h>

/*
  N-PULSES
  
  The dyadic N-pulses are defined for any {c \geq -1} and any
  non-negative degree {g \geq 2*c+1}. 
  
  For a given {g} and {c}, there are {g-c} distinct mother pulses
  {w^g_c[pix]}. They fall into two sub-classes:

  * When {pix} is in {c+1..g-c-1}, the mother spline has
    support {(0_1)}, and coincides with the Bernstein-Bézier
    polynomial {BB^g_pix}:

      { w^g_c[pix](z) = BB^g_pix(z) = choose(g,pix) z^{g-pix} (1-z)^pix }

  * When {pix} is in {0..c}, the mother spline is supported on 
    two grid cells, spanning {(0_2)}. In the first cell {(0_1)},
    it can be expressed as a combination of Bernstein-Bézier
    polynomials:
    
      { w^g_c[pix](z) = SUM { choose(k,pix)/2^{k} BB^g_{g-c+k}(z) : k = pix..c } }
    
    For {z} in {(1_2)}, the following relation holds

      { w^g_c[pix](z) = w^g_c[c-pix](2-z) }
      
    so 
      
      { w^g_c[pix](z) = SUM { choose(c-k,pix)/2^{c-k} BB^g_{g-k}(z) : k = 0..c-pix }
        
  */

unsigned int psp_pulse_mother_N_count (psp_cont_t c, psp_degree_t g);
  /* Returns the number of mother N-pulses with continuity {c}
    and degree {g}. Returns 0 if such N-pulses aren't defined or 
    aren't supported by this module. */

psp_grid_size_t psp_pulse_mother_N_supp_count
  ( psp_cont_t c, 
    psp_degree_t g, 
    psp_pulse_mother_index_t pix
  );
  /* Returns the number of unit-grid cells that comprise the support
    of the mother N-pulse with given {c,g,pix}. Fails if such N-pulses
    aren't defined or aren't supported by this module. */
    
void psp_pulse_mother_N_to_bezier
  ( psp_cont_t c, 
    psp_degree_t g, 
    psp_pulse_mother_index_t pix,
    psp_grid_pos_t x,
    double *bz
  );
  /* Stores into {bz[0..g]} the Bézier coefficients of the mother
    N-pulse with continuity {c}, degree {g}, and mother index {pix}, for
    the interval {(x _ x+1)}. Fails if such pulse is not defined or
    supported by this module. The coefficients will be non-zero only
    if {x == 0}, or if {x == 1} and {pix <= c}. */

#endif
