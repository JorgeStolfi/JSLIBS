/* Picking slicing planes.  */
/* Last edited on 2024-10-06 03:39:55 by stolfi  */

#ifndef tosl_pick_planes_H
#define tosl_pick_planes_H

#define _GNU_SOURCE
#include <stdint.h>

#include <tosl.h>
      
tosl_coord_t *tosl_pick_planes(int32_t NP, tosl_coord_t Z0, tosl_coord_t Z1, int32_t type, int32_t verbose);
  /* Returns a vector of {NP} {Z}-coordinates of planes, strictly increasing,
    starting with {Z0} and ending with {Z1}.
    
    The ends {Z0} and {Z1} must be odd, and the difference {Z1-Z0} must be at least twice {NP-1}.
    
    The {type} determines the spacing of the planes. It {type} is 0, the
    spacing will be as uniform as possible. If {type} is 1, the spacing
    will have both random and gradual variations. If {type} is 3, the
    spacing will be semicountinous in a "fractal" sense.
    
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
