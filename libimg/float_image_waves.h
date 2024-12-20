#ifndef float_image_waves_H
#define float_image_waves_H

/* Tools for generating images that are combinations of waves. */
/* Last edited on 2024-12-05 10:29:58 by stolfi */ 

#include <stdint.h>

#include <float_image.h>

double float_image_waves_eval
  ( double x,
    double y, 
    int32_t NF,
    double amp[], 
    double fx[], 
    double fy[], 
    double phase[],
    double squash
  );
  /* Computes the sum of {NF} sinusoidal waves at the point {(x,y)}.
    Each wave has amplitude {amp[k]}, frequency vectors {(fx[k],fy[k])},
    and phase {phase[k]}, for {k} in {0..NF-1}.
    
    The frequency vectors should have length less than 0.5, to avoid
    aliasing artifacts when sampled in a 1-pixel grid.
    
    If {squash} is finite and positive, the result is processed by a 
    map {F) that compresses the range {(-INF _ +INF)} to {[-1 _ +1]},
    such that {F(±squash) = ±sqrt(0.5)}. */

void float_image_waves_pick
  ( int32_t NF,
    double amp[], 
    double fx[], 
    double fy[], 
    double phase[],
    bool_t verbose
  );
  /* Chooses a set of {NF} wave functions for {float_image_test_comb_waves}.
    
    The wawes come from a series with exponentially decreasing frequency vectors 
    and amplitudes. The procedure stores the parameters in {amp[k]},
    (fx[k]}, {fy[k]}, and {phase[k]}, for {k} in {0..NF-1}.  
    
    If {verbose} is true, prints the wave parameters to {stderr}. */

#endif
