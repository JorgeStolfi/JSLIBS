#ifndef float_image_waves_H
#define float_image_waves_H

/* Tools for generating images that are combinations of waves. */
/* Last edited on 2023-01-10 20:08:03 by stolfi */ 

#define _GNU_SOURCE_
#include <stdint.h>

#include <float_image.h>

double float_image_waves_eval
  ( double x,
    double y, 
    int32_t NF,
    double amp[], 
    double fx[], 
    double fy[], 
    double phase[]
  );
  /* Computes the sum of {NF} sinusoidal waves at the point {(x,y)}.
    Each wave has amplitude {amp[k]}, frequency vectors {(fx[k],fy[k])},
    and phase {phase[k]}, for {k} in {0..NF-1}.
    
    The frequency vectors should have length less than 0.5, to avoid
    aliasing artifacts when sampled in a 1-pixel grid. */

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
