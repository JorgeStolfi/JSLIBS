/* voxb_erolate.h --- erosion and dilation of PPV arrays. */
/* Last edited on 2021-06-10 17:34:10 by jstolfi */

#ifndef voxb_erolate_H
#define voxb_erolate_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <ppv_array.h>

#define voxb_erolate_rad_MAX (20.0)
  /* Max absolute smothing radius allowed. */

void voxb_erolate_with_ball(ppv_array_t *A, double rad);
 /* Erodes or dilates the shape stored in the binary PPV array {A} by the radius {rad}.
   
   'Specifically, if {srm} is positive, erodes the shape by expanding every '0' voxel
   into a digital ball of '0' voxels.  If {srm} is negative, dilates the shape
   by expanding every '1' voxel into a ball of '1' voxels. 
   
   In either case, the ball will consist of all voxels whose centers lie
   at distance {|rad|} or less from the voxel in question (in voxels).
   This procedure runs in time proportional to {N*|rad|^3} where {N} is
   the number of voxels in {A}. */

#endif

