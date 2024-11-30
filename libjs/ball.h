#ifndef ball_H
#define ball_H

#include <stdlib.h>
#include <stdint.h>

/* Ball integrals. */
/* Last edited on 2024-11-23 06:19:27 by stolfi */
  
double ball_vol(uint32_t d);
  /* Returns the volume of the {d}-ball with radius {r}, namely
    {V(d) = V(d) = Pi^(d/2)/((d/2)!)} */

double ball_cap_vol_frac_pos(uint32_t d, double u);
  /* Returns the fraction {S(d,u)} of the volume of the unit {d}-ball
    that is contained in the slice between {x=-1} and {x=u},
    for {u} in {[-1,+1]}. */

double ball_zone_vol_frac_ang(uint32_t d, double w);
  /* Returns the fraction of the volume of the unit {d}-ball contained
    in the slice between {x=0} and {x=sin(w)}, for {w} in
    {[-PI/2,PI/2]}; namely, 
    
      {F(d,w) = (V(d-1)/V(d))*integral((cos(z))^d, z \in 0 _ w)}
      
    where {V(d)} is the volume of the unit {d}-ball. */

#endif
