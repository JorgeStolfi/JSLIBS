/* Picking slicing planes.  */
/* Last edited on 2024-10-07 06:28:07 by stolfi  */

#ifndef tosl_pick_planes_H
#define tosl_pick_planes_H

#define _GNU_SOURCE
#include <stdint.h>

#include <tosl.h>
      
tosl_coord_t *tosl_pick_planes(int32_t NP, tosl_coord_t Z0, tosl_coord_t Z1, int32_t type, int32_t verbose);
  /* Returns a vector {Zplane[0..NP-1]} of {Z}-coordinates of sicing planes, all odd and strictly increasing,
    starting with {Z0} and ending with {Z1}.
    
    The end values {Z0} and {Z1} must be odd, and the difference {Z1-Z0} must be at least twice {NP-1}.
    
    The {type} determines the spacing of the planes:
    
      Type 0: {Zplane[ip]} will be as close to affine interpolation between {Z0} and {Z1}
        as possible, namely {Z0 + ip*(Z1-Z0)/(NP-1)} except for the roundoff of
        quantization to odd integers.
        
      Type 1: the spacing will have large variations, partly random
        and partly a gradual function of {ip}.
     
      Type 2: the spacing will have large variations in a "brownian fractal" style.
      
      Type 3: the spacing will vary in sections, so that the largest spacing
        is at most four times the smallest spacing.
      
    Eacept for type 0, the value of {Zplane[ip]} will deviate significantly from
    affine interpolation between {Z0} and {Z1}.
    
    If the boolean {verbose} is true, the procedure prints various diagnostics as it 
    works. */

void tosl_pick_planes_print_stats
  ( FILE *wr,
    int32_t NP,
    tosl_coord_t Zplane[],
    tosl_coord_t Zmin,
    tosl_coord_t Zmax
  );
  /* Prints to {wr} various statistics about the planes {Z}-values {Zplane[0..NP-1]}.
    Expects that those values will be strictly increasing and contained in the
    range {Zmin..Zmax} of vertex {Z}-coordinates. */

#endif
