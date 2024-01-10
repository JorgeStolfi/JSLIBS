/* N-pulses - Compact unit-partition finite elements. */
/* Last edited on 2009-02-09 18:23:36 by stolfi */

#ifndef udg_pulse_N_H
#define udg_pulse_N_H

#include <bz_basic.h>

#include <udg_pulse.h>

/*
  N-PULSES
  
  The dyadic N-pulses are defined for any {c \geq -1} and any
  non-negative degree {g \geq 2*c+1}. 
  
  For a given {g} and {c}, there are {g-c} distinct mother pulses
  {w^g_c[pix]}. They fall into two sub-classes:

  * When {pix} is in {c+1..g-c-1}, the mother spline has
    support {(0_1)}, and coincides with the Bernstein-B�zier
    polynomial {BB^g_pix}:

      { w^g_c[pix](z) = BB^g_pix(z) = choose(g,pix) z^{g-pix} (1-z)^pix }

  * When {pix} is in {0..c}, the mother spline is supported on 
    two grid cells, spanning {(0_2)}. In the first cell {(0_1)},
    it can be expressed as a combination of Bernstein-B�zier
    polynomials:
    
      { w^g_c[pix](z) = SUM { choose(k,pix)/2^{k} BB^g_{g-c+k}(z) : k = pix..c } }
    
    For {z} in {(1_2)}, the following relation holds

      { w^g_c[pix](z) = w^g_c[c-pix](2-z) }
      
    so 
      
      { w^g_c[pix](z) = SUM { choose(c-k,pix)/2^{c-k} BB^g_{g-k}(z) : k = 0..c-pix }
        
  */

unsigned int udg_pulse_N_num_mothers (udg_cont_t c, udg_degree_t g);
  /* Returns the number of mother N-pulses with continuity {c}
    and degree {g}. Returns 0 if such N-pulses aren't defined or 
    aren't supported by this module. */

udg_grid_size_t udg_pulse_N_mother_supp_count
  ( udg_cont_t c, 
    udg_degree_t g, 
    udg_pulse_mother_index_t pix
  );
  /* Returns the number of unit-grid cells that comprise the support
    of the mother N-pulse with given {c,g,pix}. Fails if such N-pulses
    aren't defined or aren't supported by this module. */
    
void udg_pulse_N_mother_to_bezier
  ( udg_cont_t c, 
    udg_degree_t g, 
    udg_pulse_mother_index_t pix,
    udg_grid_pos_t x,
    double *bz
  );
  /* Stores into {bz[0..g]} the B�zier coefficients of the mother
    N-pulse with continuity {c}, degree {g}, and mother index {pix}, for
    the interval {(x _ x+1)}. Fails if such pulse is not defined or
    supported by this module. The coefficients will be non-zero only
    if {x == 0}, or if {x == 1} and {pix <= c}. */

#endif
