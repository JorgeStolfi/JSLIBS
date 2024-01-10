/* B-pulses - B-Spline-like. */
/* Last edited on 2009-02-09 18:22:32 by stolfi */

#ifndef udg_pulse_B_H
#define udg_pulse_B_H

#include <bz_basic.h>

#include <udg_pulse.h>

/*
  B-PULSES
  
  The dyadic B-pulses are defined for any {c \geq -1}. 
  
  For a given {c}, there is exactly one mother B-pulse, with degree
  {g = c+1} and index {pix = 0}. It is supported on {msz = c+2} intervals, 
  spanning the interval {(0 _ msz)}.
  
  For {c = -1} the mother B-pulse is the unit constant function
  in the interval {(0_1)}.  For {c = 0} the B-pulse is the 
  same as the B-pulse and the H-pulse (a triangular function
  on {(0_2)}. For {c \geq 2} the B-pulses are the B-spline 
  basis. */

unsigned int udg_pulse_B_num_mothers (udg_cont_t c, udg_degree_t g);
  /* Returns the number of mothers for B-pulses with continuity {c}
    and degree {g}: namely, 1 if {g = c+1}, 0 otherwise. */

udg_grid_size_t udg_pulse_B_mother_supp_count
  ( udg_cont_t c, 
    udg_degree_t g, 
    udg_pulse_mother_index_t pix
  );
  /* Returns the support count of the mother B-pulse with given
    {c,g,pix}: namely, {g+1 = c+2} if {g == c+1}, undefined
    otherwise. */
    
void udg_pulse_B_mother_to_bezier
  ( udg_cont_t c, 
    udg_degree_t g, 
    udg_pulse_mother_index_t pix,
    udg_grid_pos_t x,
    double *bz
  );
  /* Stores into {bz[0..g]} the Bézier coefficients of the mother
    B-pulse with continuity {c}, degree {g}, and mother {pix}, for the
    interval {(x _ x+1)}. Fails if such pulse is not defined or
    supported by this module. The coefficients will be non-zero only
    if {x == 0}, or if {x == -1} and {pix <= c}. */

#endif
